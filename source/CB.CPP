#include "CB.h"

static const double epsV = 1e-6;		// Погрешность определения скорости,м/с
static const double epsM = 1e-6;		// Погрешность определения массы,кг
static const double epsT = 1e-8;		// Погрешность определения времени,с
static const double epsR = 1e-5;		// Погрешность определения радиус-вектора,м
static const double epsTetta = 1e-7;	// Погрешность определения угла наклона траектории,рад
static const double eps_extr = 1e-2;	// Точность решения параметрической оптимизации

double CB::tetta_control(double t)
{
	if (t <= m_tv)			// Вертикальный участок полёта
		return M_PI_2;
	else if (t <= m_t1)		// 1-ый АУТ (тангаж - линейная функция)
		return M_PI_2 + m_tetta_point * (t - m_tv);	
	else if (t < m_t1 + m_tPP)
		return 0;			// Пассивный участок полёта - нет тяги
	else
		return m_tetta2;	// 2-ой АУТ (тангаж - постоянная функция)
}

double CB::get_Tetta(double params[])
{
	if (abs(get_R(params[2], params[3]) * get_V(params[0], params[1])) < 0.00001)
		return M_PI_2;
	else
		return asin((params[2] * params[0] + (params[3] + R_M) * params[1]) / 
			(get_R(params[2], params[3]) * get_V(params[0], params[1])));
}

void CB::Step(double & t, double dt, double params[], double old_params[])
{
	double K[5][4] = { 0 };			// Параметры К (производные) в методе РК4 (собраны в матрицу)
	double params_temp[5] = { 0 };	// Копия интегрируемых параметров (невмешательство в основные)
	// Снимаем копию с основных параметров - расчёт будет для этой копии
	for (int i = 0; i < 5; ++i)
		params_temp[i] = params[i];
	double delta_t[4] = { 0, dt / 2 , dt / 2 , dt };	// Подшаги в методе РК4
	// Непосредственное интегрирование - расчёт коэф-тов К в методе РК4
	for (int i = 0; i <= 3; ++i)
	{
		K[0][i] = dVx_dt(params_temp[4], t + delta_t[i], params_temp[2], params_temp[3]);
		K[1][i] = dVy_dt(params_temp[4], t + delta_t[i], params_temp[2], params_temp[3]);
		K[2][i] = dx_dt(params_temp[0]);
		K[3][i] = dy_dt(params_temp[1]);
		K[4][i] = dm_dt(m_Betta);

		if (i != 3)		// При i = 3 будет выход за размеры массива delta_t
		{
			params_temp[0] = params[0] + K[0][i] * delta_t[i + 1];
			params_temp[1] = params[1] + K[1][i] * delta_t[i + 1];
			params_temp[2] = params[2] + K[2][i] * delta_t[i + 1];
			params_temp[3] = params[3] + K[3][i] * delta_t[i + 1];
			params_temp[4] = params[4] + K[4][i] * delta_t[i + 1];
		}
	}
	// Запоминаем параметры до интегрирования
	for (int i = 0; i < 5; ++i)
		old_params[i] = params[i];
	// Расчёт параметров после интегрирования методом РК4 на один шаг
	for (int i = 0; i < 5; ++i)
		params[i] += (K[i][0] + 2 * K[i][1] + 2 * K[i][2] + K[i][3]) * dt / 6;
	// Увеличение времени на величину шага
	t += dt;
}

CB::CB::Nevyazki CB::RungeKutt(Hkon h_index, Pkon P_index)
{
	double t = 0;			// Время от начала интегрирования,с
	double dt = 0.1;		// Шаг интегрирования,с
	double t_old = 0;		// Время на предыдущем шаге интегрирования,с
	double t_prev[11]{ 0 }, t_prevT1[11]{ 0 }, t_prevT2[11]{ 0 };// Массивы для времён дробления шага интегрирования,с
	double delta_V;			// Невязка по скорости,м/с
	double params[5] = { 0 };		// Параметры интегрирования (см. ниже)
	double old_params[5] = { 0 };	// Параметры интегрирования с предыдущего шага
	m_P = arr_P[P_index], m_Betta = arr_Betta[P_index];	// Выбор тяги и секундного расхода
	m_h = arr_h[h_index][P_index];						// Выбор высоты выведения
	double R_circle = R_M + m_h;						// Радиус круговой орбиты выведения,м
	double V_orbit = sqrt(mu_M / R_circle);				// Скорость КА на круговой орбите,м/с
	int count = 0, countT1 = 0, countT2 = 0;			// Счётчики для дробления шага
	int PrintStep = 0;									// Номер шага печати (номер строки вывода)
	int Print_Time = 1;									// Шаг печати,с
	
	bool isDiv;		// Логическая переменная - дробить ли шаг?

	params[0] = m_Vx0;		// Vx,м/с;
	params[1] = m_Vy0;		// Vy,м/с;
	params[2] = m_x0;		// x,м;
	params[3] = m_y0;		// y,м;
	params[4] = m_m0;		// m,кг;

	if (isPrint)
	{
		PrintReport(true, false, false, t, params);	// Печать шапки таблицы
		PrintAllParams(true, false, false, t, params);
	}
	// Цикл интегрирования траектории (от старта СВ до выполнения условия останова)
	while (true)
	{
		isDiv = true;
		// Запись в файл
		if ((abs(t - Print_Time * PrintStep) < 1e-8) && (isPrint))
		{
			PrintReport(false, false, false, t, params);	// Печать параметров (1 строка для текущего t)
			PrintAllParams(false, false, false, t, params);
			++PrintStep;	// Переход на след. шаг печати (чтобы не печатать два раза одну строку при возврате назад)
		}
		t_old = t;		// Запоминается текущее время - до очередного шага интегрирования
		Step(t, dt, params, old_params);		// Один шаг методом интегрирования (РК4)
		delta_V = get_deltaV(get_V(params[0], params[1]), V_orbit);	// Разница текущей скорости и требуемой (орбитальной)
		// Условие останова (набор СВ орбитальной скорости с заданной точностью)
		if (fabs(delta_V) < epsV)	// Достижение заданной (epsV) точности
		{
			break;		// Прерывание цикла while (окончание интегрирования)
		}
		if (delta_V * get_deltaV(get_V(old_params[0], old_params[1]), V_orbit) < 0)
		{
			// Возврат на предыдущий шаг (чтобы затем подойти с меньшим шагом)
			for (int i = 0; i < 5; i++)
				params[i] = old_params[i];
			t -= dt;	// В т.ч. необходимо уменьшить время на величину шага
			dt /= 10;	// Дробление шага (деление в 10 раз)
			isDiv = false;	// Не должно быть одновременного дробления и по скорости и по массе!!!
		}
		// Дробление шага при израсходовании топлива (точно определяем момент выключения двигателя)
		StepFragment(get_deltaM(params[4]), get_deltaM(old_params[4]), epsM, isDiv, params, old_params, t, dt, t_prev, count, true, P_index);
		// Дробление шага при достижении t1 (точно определяем момент выключения двигателя)
		StepFragment(t - m_t1, t_old - m_t1, epsT, isDiv, params, old_params, t, dt, t_prevT1, countT1, true, P_index);
		// Дробление шага при достижении t1 + tPP (точно определяем момент включения двигателя)
		StepFragment(t - (m_t1 + m_tPP), t_old - (m_t1 + m_tPP), epsT, isDiv, params, old_params, t, dt, t_prevT2, countT2, false, P_index);
	}
	if (isPrint)	// Запись в файл параметров в конце траектории
	{
		PrintReport(false, true, true, t, params);
		PrintAllParams(false, true, true, t, params);
	}
	// Запоминаем пар-ры в конце траектории
	m_Vx = params[0], m_Vy = params[1], m_x = params[2], m_y = params[3],
		m_mk = params[4], m_tk = t;
	// Создаём объект структуры для записи невязок
	Nevyazki discrepancies;
	discrepancies.V = get_deltaV(get_V(m_Vx, m_Vy), V_orbit);		// Невязка по скорости
	discrepancies.r = get_deltaR(get_R(m_x, m_y), R_circle);		// Невязка по радиус вектору
	discrepancies.Tetta = get_deltaTetta(get_Tetta(params));		// Невязка по углу наклона траектории
	discrepancies.Massa = get_deltaM(params[4]);					// Невязка по массе (остатки топлива)
	return discrepancies;
}

void CB::StepFragment(double current, double old, double eps, bool isDiv, double params[], double old_params[], 
	double & t, double & dt, double t_prev[], int & count, bool CutOff, int P_index)
{
	if ((current * old < 0) && (fabs(current) > eps) && (isDiv))	// Исследуемый параметр переходит через ноль
	{
		for (int i = 0; i < 5; i++)
			params[i] = old_params[i];	// Возврат к параметрам на предыдущем шаге
		t -= dt;			// Возврат к времени на пердыдущем шаге
		t_prev[count] = t;	// Запоминается время интегрирования, на котором происходило дробление шага
		++count;			// Счётчик отсчитывает, сколько раз было проведено дробление шага
		dt /= 10;			// Шаг уменьшается в 10 раз
	}
	// Достигли нужной точности определения массы
	if ((current * old < 0) && (fabs(current) < eps))	// Исследуемый пар-р переходит через 0 и дости. заданная точность
	{
		if (CutOff)		// Если отключение двигателя
		{
			m_P = 0;		// Тяга двигателя теперь будет только нулевая
			m_Betta = 0;	// Как и секундный массовый расход - топлива больше нет
		} else {		// Включение двигателя
			m_P = arr_P[P_index];		
			m_Betta = arr_Betta[P_index];	// Выбор тяги и секундного расхода
		}
		if (isPrint)	// Запись в файл момента выключения/включения ДУ
		{
			PrintReport(false, true, false, t, params);
			PrintAllParams(false, true, false, t, params);
		}
	}
	// Возвращение к исходной величине шага интегрирования - осуществляется постепенно (как и дробление)
	if ((dt < 0.1) && (fabs(t - (t_prev[count - 1] + dt * 10)) < 1e-9) && (isDiv))
	{
		for (int i = count; i > 0; --i)	// В счётчике кол-во дроблений - будем уменьшать счётчик, восстанавливая шаг
		{
			if (fabs(t - (t_prev[count - i] + dt * pow(10, i))) < eps)
			{
				dt *= pow(10, i);	// Увеличение шага в 10 раз
				count -= i;			// Уменьшение счётчика
				break;				// 2 раза на одном шаге интегрирования увеличения быть не должно
			}
		}
	}
}

void CB::InverseMatrix(double A[2][2], double A_1[2][2])
{
	double det_A = 0;	// Определитель исходной матрицы
	det_A = A[0][0] * A[1][1] - A[1][0] * A[0][1];	// Расчёт определителя
	if (det_A == 0)		// Если определитель = 0 -> обратной матрицы не существует
	{
		std::cout << "ERROR: det_A = 0" << std::endl;
	} else {
		double M[2][2]{ 0 };		// Матрица алгебраических дополнений (Матрица миноров с изм. знаками на второст. диагонали)
		double M_T[2][2]{ 0 };		// Транспонированная матрица алгебраических дополнений
		M[0][0] = A[1][1], M[0][1] = -A[1][0];
		M[1][0] = -A[0][1], M[1][1] = A[0][0];
		M_T[0][0] = M[0][0], M_T[0][1] = M[1][0];
		M_T[1][0] = M[0][1], M_T[1][1] = M[1][1];
		// Рассчитываем все элементы обратной матрицы (двумя циклами - по строкам и столбцам)
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j)
				A_1[i][j] = M_T[i][j] / det_A;	// Элемент обратной матрицы
	}
}

void CB::SolveKZ()
{
	Nevyazki nevyazki, nevyazki_temp;		// Результаты работы метода интегрирования - невязки в конце траектории
	double d_tetta_point = 0, d_tetta2 = 0;	// Приращения для параметров управления (тетта с точкой и t1)
	double J[2][2]{ 0 }, J_1[2][2]{ 0 };	// Матрица Якоби и её обратная матрица (пока все элементы нулевые)
	double dX[2]{ 0 }, dY[2]{ 0 };			// Приращения на очередной итерации решения краевой задачи и невязки
	double R_k, R_0, Tetta_k, Tetta_0;		// Радиус и угол наклона траекториии на "эталонной" траектории и траектории при наличии приращения

	// Расчёт дельт - приращений для расчёта частных производных матрицы Якоби

	d_tetta_point = m_tetta_point / 1e+5;
	d_tetta2 = m_tetta2 / 1e+5;

	while (true)
	{
		// Первое интегрирование (без приращений параметров t_point и t1)

		nevyazki = RungeKutt(h_High, P_High);

		R_0 = get_R(m_x, m_y);		// Значение радиус-вектора в конце траектории
		Tetta_0 = nevyazki.Tetta;	// Значение угла наклона траектории --//--
		dY[0] = nevyazki.r, dY[1] = nevyazki.Tetta;	// Вектор невязок (радиус-вектор и Тетта)

		if ((fabs(dY[0]) < epsR) && (fabs(dY[1]) < epsTetta))	// Условие прерывания бесконечного цикла
		{
			//isPrint = true;		// Запись параметров траектории только при выполнении всех ГУ
			nevyazki = RungeKutt(h_High, P_High);	// "Чистовое" интегрирование
			m_dm = nevyazki.Massa;					// Остатки топлива при чистовом расчёте
			break;	// Выход из цикла
		}

		// Второе интегрирование и определение невязок

		m_tetta_point += d_tetta_point;		// Прибавляем малое приращение к tetta_point
		nevyazki_temp = RungeKutt(h_High, P_High);	// Интегрируем траекторию с новым значением tetta_point
		R_k = get_R(m_x, m_y);			// Значение радиус-вектора в конце траектории
		Tetta_k = nevyazki_temp.Tetta;	// Значение угла наклона траектории --//--

		// Первый СТОЛБЕЦ матрицы Якоби

		J[0][0] = (R_k - R_0) / d_tetta_point;	// Первый (a11) элемент матрицы Якоби
		J[1][0] = (Tetta_k - Tetta_0) / d_tetta_point;	// Третий (a21) элемент матрицы Якоби
		m_tetta_point -= d_tetta_point;		// Убираем малое приращение от tetta_point (возврат к исходному значению)

		// Третье интегрирование и определение невязок

		m_tetta2 += d_tetta2;				// Прибавляем малое приращение к tetta2
		nevyazki_temp = RungeKutt(h_High, P_High); // Интегрируем траекторию с новым значением tetta2
		R_k = get_R(m_x, m_y);				// Значение радиус-вектора в конце траектории
		Tetta_k = nevyazki_temp.Tetta;		// Значение угла наклона траектории --//--

		// Второй СТОЛБЕЦ матрицы Якоби

		J[0][1] = (R_k - R_0) / d_tetta2;			// Второй (a12) элемент матрицы Якоби
		J[1][1] = (Tetta_k - Tetta_0) / d_tetta2;	// Четвёртый (a22) элемент матрицы Якоби
		m_tetta2 -= d_tetta2;						// Убираем малое приращение от tetta2 (возврат к исходному значению)

		InverseMatrix(J, J_1);		// Обращение матрицы Якоби
		dX[0] = J_1[0][0] * dY[0] + J_1[0][1] * dY[1];		// Приращение для tetta_point
		dX[1] = J_1[1][0] * dY[0] + J_1[1][1] * dY[1];		// Приращение для tetta2

		// Определение новых значений параметров управления

		m_tetta_point -= dX[0] * 1.0;		// Коэффициент вводится для более точного (аккуратного) поиска
		m_tetta2 -= dX[1] * 1.0;			// Чем коэффициент меньше - тем больше итераций решения КЗ
	}
}

void CB::Optimization()
{
	double d_t1 = 0, d_tPP = 0;							// Приращения для параметров управления (t1 и t_ПУТ)
	double alpha = 1e-3;								// Числовой коэффициент (задаётся его начальное приближение)
	const double C1 = 1.2, C2 = 1 / C1;					// Числовые коэффициенты, подчиняющиеся условиям: C1 [0,1], C2 > 1
	const double I[2][2]{ {1, 0}, {0, 1} };				// Единичная матрица 2х2
	double H[2][2]{ 0 }, grad[2]{ 1e+4, 1e+4 };			// Матрица Гессе и градиент
	double H1[2][2]{ 0 }, H_1[2][2]{ 0 };				// Матрица для определения след. приближ. и её обратная (она и нужна)
	double g_i[2]{ 0 };									// Минимизируемая функция (израсходованное топливо за время полёта)
	double m_k, m_0, m_k2[2];							// Значения массы в конце траектории для расчёта частных производных
	double dm_dtPP[2], dm_dt1[2];						// Первые частные производные - для расчёта вторых ЧП
	double df[2][2]{ 0 };
	int i = 0;											// Счётчик для печати итераций в консоль
	//Расчёт дельт - приращений для расчёта частных производных матрицы Гессе
	//d_t1 = m_t1 / 1e+5;
	d_t1 = 0.01;
	//d_tPP = m_tPP / 1e+5;
	d_tPP = 0.01;
	// Поиск оптимального решения
	while (true)
	{
		// Первое решение КЗ (без малых приращений параметров t_ПУТ и tetta2)
		SolveKZ();										// Решение краевой задачи
		m_0 = m_mk;										// Конечная масса СВ,кг
		g_i[0] = m_m0 - m_mk;							// Затраты топлива,кг (минимизируемая функция)
// РАСЧЁТ ВЕКТОРА ГРАДИЕНТА
		// Определение производной dm/dt_ПУТ
		
		grad[0] = dg_dtPP(d_tPP / 2);					// Первый элемент вектора градиента
		
		// Определение производной dm/dt1
		
		grad[1] = dg_dt1(d_t1 / 2);						// Второй элемент вектора градиента
		
		// Конец расчёта вектора градиента
		std::cout << i << "\tGrad_1 = " << grad[0] << "\tGrad_2 = " << grad[1] << "\tПотрачено = " << g_i[0] << std::endl;
		i++;
// ПРОВЕРКА УСЛОВИЯ ОСТАНОВА ПОИСКА ОПТИМАЛЬНОГО РЕШЕНИЯ
		if (sqrt(pow(grad[0], 2) + pow(grad[1], 2)) < eps_extr)	// Условие прерывания бесконечного цикла
		{
			isPrint = true;		// Запись параметров траектории только при выполнении всех ГУ
			SolveKZ();			// Решение краевой задачи (чистовое - запись файлов)
			break;				// Выход из цикла
		}
// РАСЧЁТ МАТРИЦЫ ГЕССЕ
		// Определение второй частной производной d2m/dt_ПУТ2
		
		H[0][0] = d2g_dtPP2(d_tPP);

		// Определение второй частной производной d2m/dt_ПУТ*dtetta2

		H[0][1] = d2g_dtPPt1(d_tPP);

		// Определение второй частной производной d2m/dtetta2*dt_ПУТ
		
		//H[1][0] = d2g_dt1tPP(d_t1);
		H[1][0] = H[0][1];
		
		// Определение второй частной производной d2m/dt1

		H[1][1] = d2g_dt12(d_t1);

		// Конец расчёта матрицы Гессе
// Определение следующего приближения к решению
		do
		{
			// Расчёт матрицы след. прибл. в соответствии с алг. Левенберга-Марквардта
			for (int i = 0; i < 2; ++i)
				for (int j = 0; j < 2; ++j)
					H1[i][j] = H[i][j] + alpha * I[i][j];
			// Обращение вышеуказанной матрицы
			InverseMatrix(H1, H_1);
			// Новые значения параметров управления
			m_tPP -= (H_1[0][0] * grad[0] + H_1[0][1] * grad[1]);
			m_t1 -= (H_1[1][0] * grad[0] + H_1[1][1] * grad[1]);
			// Решение краевой задачи
			SolveKZ();
			// Затраты топлива,кг (минимизируемая функция)
			g_i[1] = m_m0 - m_mk;
			if (g_i[1] >= g_i[0])
				alpha *= C2;		// Изменение коэффициента alpha в соотв. с алг. Левенберга-Марквардта
		} while (g_i[1] >= g_i[0]);	// Цикл повтор-ся пока затраты не станут меньше, чем на пред. итерации
		alpha *= C1;	// Изменение коэффициента alpha в соотв. с алг. Левенберга-Марквардта
	}
}

double CB::dg_dt1(double dt)
{
	double dm_left, dm_right;
	m_t1 -= dt;					// Малое приращение параметра влево
	SolveKZ();					// Решение краевой задачи
	dm_left = m_m0 - m_mk;		// Разница масс в конце при малом прир. влево
	m_t1 += 2 * dt;				// Малое приращение параметра вправо (+1 до "нуля" и ещё +1 вправо)
	SolveKZ();					// Решение краевой задачи
	dm_right = m_m0 - m_mk;		// Разница масс в конце при малом прир. вправо
	m_t1 -= dt;					// Возвращение параметра на "ноль"
	return (dm_right - dm_left) / (2 * dt);		// Расчёт ЧП методом центральной разности
}

double CB::dg_dtPP(double dt)
{
	double dm_left, dm_right;
	m_tPP -= dt;
	SolveKZ();
	dm_left = m_m0 - m_mk;
	m_tPP += 2 * dt;
	SolveKZ();
	dm_right = m_m0 - m_mk;
	m_tPP -= dt;
	return (dm_right - dm_left) / (2 * dt);
}

double CB::d2g_dt12(double dt)
{
	double df[2][2]{ 0 }, dm_dt1[2]{ 0 };
	m_t1 -= 2 * dt;
	SolveKZ();
	df[0][0] = m_m0 - m_mk;
	m_t1 += 2 * dt;
	SolveKZ();
	df[0][1] = m_m0 - m_mk;
	dm_dt1[0] = (df[0][1] - df[0][0]) / (2 * dt);
	m_t1 += 2 * dt;
	SolveKZ();
	df[1][0] = m_m0 - m_mk;
	m_t1 -= 2 * dt;
	//SolveKZ();
	df[1][1] = df[0][1];
	dm_dt1[1] = (df[1][0] - df[1][1]) / (2 * dt);
	return (dm_dt1[1] - dm_dt1[0]) / (2 * dt);
}

double CB::d2g_dtPP2(double dt)
{
	double df[2][2]{ 0 }, dm_dtPP[2]{ 0 };
	m_tPP -= 2 * dt;
	SolveKZ();
	df[0][0] = m_m0 - m_mk;
	m_tPP += 2 * dt;
	SolveKZ();
	df[0][1] = m_m0 - m_mk;
	dm_dtPP[0] = (df[0][1] - df[0][0]) / (2 * dt);
	m_tPP += 2 * dt;
	SolveKZ();
	df[1][0] = m_m0 - m_mk;
	m_tPP -= 2 * dt;
	//SolveKZ();
	df[1][1] = df[0][1];
	dm_dtPP[1] = (df[1][0] - df[1][1]) / (2 * dt);
	return (dm_dtPP[1] - dm_dtPP[0]) / (2 * dt);
}

double CB::d2g_dt1tPP(double dt)
{
	double df[2][2]{ 0 }, dm_dtPP[2]{ 0 };
	m_t1 -= dt;
	m_tPP -= dt;
	SolveKZ();
	df[0][0] = m_m0 - m_mk;
	m_tPP += 2 * dt;
	SolveKZ();
	df[0][1] = m_m0 - m_mk;
	m_tPP -= dt;
	dm_dtPP[0] = (df[0][1] - df[0][0]) / (2 * dt);
	m_t1 += 2 * dt;
	m_tPP -= dt;
	SolveKZ();
	df[1][0] = m_m0 - m_mk;
	m_tPP += 2 * dt;
	SolveKZ();
	df[1][1] = m_m0 - m_mk;
	m_tPP -= dt;
	dm_dtPP[1] = (df[1][1] - df[1][0]) / (2 * dt);
	m_t1 -= dt;
	return (dm_dtPP[1] - dm_dtPP[0]) / (2 * dt);
}

double CB::d2g_dtPPt1(double dt)
{
	double df[2][2]{ 0 }, dm_dt1[2]{ 0 };
	m_tPP -= dt;
	m_t1 -= dt;
	SolveKZ();
	df[0][0] = m_m0 - m_mk;
	m_t1 += 2 * dt;
	SolveKZ();
	df[0][1] = m_m0 - m_mk;
	m_t1 -= dt;
	dm_dt1[0] = (df[0][1] - df[0][0]) / (2 * dt);
	m_tPP += 2 * dt;
	m_t1 -= dt;
	SolveKZ();
	df[1][0] = m_m0 - m_mk;
	m_t1 += 2 * dt;
	SolveKZ();
	df[1][1] = m_m0 - m_mk;
	m_t1 -= dt;
	dm_dt1[1] = (df[1][1] - df[1][0]) / (2 * dt);
	m_tPP -= dt;
	return (dm_dt1[1] - dm_dt1[0]) / (2 * dt);
}

void CB::PrintReport(bool isFirstEnter, bool isExtr, bool isLast, double t, double params[])
{
	int ind = 20;			// Отступ в текстовых файлах
	if (isFirstEnter)		// При первом вызове метода - открывается файл и записывается шапка
	{
		fout.open("ReportTable.txt");
		fout.clear();
		fout.setf(std::ios_base::fixed);
		fout.width(40);
		fout << " Результаты расчёта траектории выведения на круговую орбиту ИСЛ с высотой h = " <<
			int(m_h / 1000) << " км (при P = " << int(m_P) << " Н)" << std::endl;
		fout << std::setw(ind) << "t,c" << std::setw(ind) << "m,кг" << std::setw(ind) << "Vx,м/c" <<
			std::setw(ind) << "Vy,м/c" << std::setw(ind) << "x,км" << std::setw(ind) << "y,км" <<
			std::setw(ind) << "h,км" << std::setw(ind) << "V,м/c" << std::setw(ind) << "r,км" <<
			std::setw(ind) << "tetta,град" << std::setw(ind) << "Tetta_c,град" << std::setw(ind) <<
			"alpha,град" << std::setw(ind) << "fi,град" << std::endl;
	} else {		// При дальнейших вызовах происходит запись параметров
		if (isExtr)
		{
			fout << std::setw(ind) << std::setprecision(6) << t << std::setw(ind) << std::setprecision(1) << params[4] << std::setw(ind) <<
				std::setprecision(4) << params[0] << std::setw(ind) << params[1] << std::setw(ind) << params[2] / 1e+3 << std::setw(ind) <<
				params[3] / 1e+3 << std::setw(ind) << get_h(params[2], params[3]) / 1e+3 << std::setw(ind) << get_V(params[0], params[1]) <<
				std::setw(ind) << get_R(params[2], params[3]) / 1e+3 << std::setw(ind) << RadianToDegree(tetta_control(t)) << std::setw(ind) <<
				RadianToDegree(get_Tetta_C(params)) << std::setw(ind) << RadianToDegree(get_alpha(t, params)) << std::setw(ind) <<
				RadianToDegree(get_fi(params)) << std::endl;
		} else {
			fout << std::setw(ind) << std::setprecision(0) << t << std::setw(ind) << std::setprecision(1) << params[4] << std::setw(ind) <<
				std::setprecision(2) << params[0] << std::setw(ind) << params[1] << std::setw(ind) << params[2] / 1e+3 << std::setw(ind) <<
				params[3] / 1e+3 << std::setw(ind) << get_h(params[2], params[3]) / 1e+3 << std::setw(ind) << get_V(params[0], params[1]) <<
				std::setw(ind) << get_R(params[2], params[3]) / 1e+3 << std::setw(ind) << RadianToDegree(tetta_control(t)) << std::setw(ind) <<
				RadianToDegree(get_Tetta_C(params)) << std::setw(ind) << RadianToDegree(get_alpha(t, params)) << std::setw(ind) <<
				RadianToDegree(get_fi(params)) << std::endl;
		}
	}
	if (isLast)
		fout.close();
}

void CB::PrintAllParams(bool isFirstEnter, bool isExtr, bool isLast, double t, double params[])
{
	int ind = 20;			// Отступ в текстовых файлах
	if (isFirstEnter)		// При первом вызове метода - открывается файл и записывается шапка
	{
		fout2.open("ReportAllParams.txt");
		fout2.clear();
		fout2.setf(std::ios_base::fixed);
		fout2.width(40);
		fout2 << " Результаты расчёта траектории выведения на круговую орбиту ИСЛ с высотой h = " <<
			int(m_h / 1000) << " км (при P = " << int(m_P) << " Н)" << std::endl;
		fout2 << "\tt1 = " << m_t1 << " c" << "\ttetta_point = " << m_tetta_point << " рад/c" << "\tt_ПУТ = " << m_tPP <<
			" c" << "\ttetta2 = " << m_tetta2 << " рад (" << RadianToDegree(m_tetta2) << " град)" << std::endl;
		fout2 << std::setw(ind) << "t,c" << std::setw(ind) << "m,кг" << std::setw(ind) << "Vx,м/c" <<
			std::setw(ind) << "Vy,м/c" << std::setw(ind) << "x,км" << std::setw(ind) << "y,км" <<
			std::setw(ind) << "h,км" << std::setw(ind) << "V,м/c" << std::setw(ind) << "r,км" <<
			std::setw(ind) << "tetta,град" << std::setw(ind) << "Tetta,град" << std::setw(ind) << "Tetta_c,град" <<
			std::setw(ind) << "alpha,град" << std::setw(ind) << "fi,град" << std::endl;
	}
	else {		// При дальнейших вызовах происходит запись параметров
		fout2 << std::setw(ind) << t << std::setw(ind) << params[4] << std::setw(ind) << params[0] <<
			std::setw(ind) << params[1] << std::setw(ind) << params[2] / 1e+3 << std::setw(ind) << params[3] / 1e+3 <<
			std::setw(ind) << get_h(params[2], params[3]) / 1e+3 << std::setw(ind) << get_V(params[0], params[1]) << std::setw(ind) <<
			get_R(params[2], params[3]) / 1e+3 << std::setw(ind) << RadianToDegree(tetta_control(t)) << std::setw(ind) <<
			RadianToDegree(get_Tetta(params)) << std::setw(ind) << RadianToDegree(get_Tetta_C(params)) << std::setw(ind) <<
			RadianToDegree(get_alpha(t, params)) << std::setw(ind) << RadianToDegree(get_fi(params)) << std::endl;
	}
	if (isLast)
		fout2.close();
}

void CB::Calculate()
{
	double dm, dm_left, dm_right, dmdmconst, dmdWist;
	// Установка первого приближения
	m_t1 = 410.0;							// Время окончания первого АУТ,с
	m_tetta_point = -0.0045;				// Производная угла тангажа на 1-м АУТ,рад/с
	m_tPP = 800;							// Продолжительность ПУТ,с
	m_tetta2 = DegreeToRadian(-35);			// Значение угла тангажа на 2-м АУТ,рад

	Optimization();							// Решение задачи параметрической оптимизации
	dm = m_dm;								// Остатки топлива при номин.парам-х
	
	// Вычисление производной dm_pg/dm_const (центральная разность)

	m_mconst -= 1;							// Уменьшение массы конструкции на 1 кг
	m_m0 -= 1;								// Уменьшение общей начальной массы на 1 кг
	Optimization();							// Решение задачи параметрической оптимизации
	dm_left = m_dm;							// Остатки топлива при приращении влево

	m_mconst += 2;							// Увеличение массы конструкции на 1 кг относит. исходной
	m_m0 += 2;								// Увеличение общей начальной массы на 1 кг относит. исходной
	Optimization();							// Решение задачи параметрической оптимизации
	dm_right = m_dm;						// Остатки топлива при приращении вправо
	m_mconst -= 1;							// Возврат на "ноль" массы конструкции
	m_m0 -= 1;								// Возврат на "ноль" общей начальной массы
	
	dmdmconst = (dm_right - dm_left) / 2;	// Вычисление производной (центральная разность)
	std::cout << "dm_pg/dm_const = " << dmdmconst << std::endl;	// Вывод на экран

	// Вычисление производной dm_pg/dW_ist (центральная разность)

	W_ist -= 1;								// Уменьшение скорости истечения на 1 м/с
	// Перерасчёт секундного расхода топлива с новой скоростью истечения
	arr_Betta[0] = arr_P[0] / W_ist;
	arr_Betta[1] = arr_P[1] / W_ist;
	Optimization();							// Решение задачи параметрической оптимизации
	dm_left = m_dm;							// Остатки топлива при приращении влево

	W_ist += 2;								// Увеличение скорости истечения на 1 м/с относит. исходной
	// Перерасчёт секундного расхода топлива с новой скоростью истечения
	arr_Betta[0] = arr_P[0] / W_ist;
	arr_Betta[1] = arr_P[1] / W_ist;
	Optimization();							// Решение задачи параметрической оптимизации
	dm_right = m_dm;						// Остатки топлива при приращении вправо
	W_ist -= 1;								// Возврат на "ноль" скорости истечения
	// Перерасчёт секундного расхода топлива к исходному состоянию
	arr_Betta[0] = arr_P[0] / W_ist;
	arr_Betta[1] = arr_P[1] / W_ist;
	dmdWist = (dm_right - dm_left) / 2;		// Вычисление производной (центральная разность)
	std::cout << "dm_pg/dWist = " << dmdWist << std::endl;	// Вывод на экран
}