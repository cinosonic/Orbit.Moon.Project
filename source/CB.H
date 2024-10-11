#pragma once						// Предотвращение повторного подключения одинаковых библиотек
#define _USE_MATH_DEFINES			// Использование математических констант
#include <iomanip>					// Библиотека для форматирования вывода
#include <iostream>					// Библиотека для работы с потоками ввода-вывода
#include <fstream>					// Библиотека с интерфейсом для чтения/записи из/в файл
#include <math.h>					// Математическая библиотека

class CB
{
private:
	double W_ist = 3490;								// Скорость истеч. проду. сгор. из сопла двиг.,м/с
	double m_P, m_Betta, m_h;		// Тяга, секундный расход и высота (инициализируются в функции RungeKutt)
	const double m_x0 = 0;								// Координата х в начальный момент времени t_0,м
	const double m_y0 = 0;								// Координата y в начальный момент времени t_0,м
	const double m_Vx0 = 0;								// Проекция скорости Vх в нач. мом. времени t_0,м/с
	const double m_Vy0 = 0;								// Проекция скорости Vy в нач. мом. времени t_0,м/с
	double m_m0 = 3015;									// Начальная масса СВ,кг
	const double m_mT_max = 1550;						// Предельная масса рабоч. запаса топлива ДУ СВ,кг
	double m_mconst = 615;								// Масса конструкции СВ,кг
	const double m_tv = 14;								// Продолжительность вертикального участка выведения,с

	const double mu_M = 4.903e+12;						// Гравитационная постоянная Луны,м^3/с^2
	const double R_M = 1738e+3;							// Радиус Луны,м
	// Vx, Vy, x, y, m в конечный момент времени t_k
	double m_Vx = 0, m_Vy = 0, m_x = 0, m_y = 0,
		m_mk = 0, m_tk = 0;
	// Вектор управления (кусочно-линейный с двумя активными участками)
	double m_tetta2 = 0, m_tetta_point = 0;				// Знач. угла танг. на 2-м АУТ (рад) и производная на 1-м (рад/с)
	double m_t1 = 0, m_tPP = 0;							// Время первого активного участка и время ПУТ,c
	double m_dm = 0;									// Масса остатков топлива,кг
	bool isPrint = false;								// Логическая переменная для записи в файл нужной траектории

	std::fstream fout, fout2;							// Переменные для записи в файлы
	enum Hkon
	{
		h_Low = 0,		// Индекс меньшей высоты
		h_High = 1,		// Индекс большей высоты
	}; 
	enum Pkon
	{
		P_Low = 0,		// Индекс меньшей тяги (и секундного расхода)
		P_High = 1,		// Индекс большей тяги (и секундного расхода)
	};
	const double arr_P[2] = { 8500, 10650 };				// Тяга СВ (j = 1, 2),Н
	double arr_Betta[2] = { arr_P[0] / W_ist, arr_P[1] / W_ist };	// Секундный расход топлива (j = 1, 2),кг/с
	const double arr_h[2][2] = { {125000, 195000}, {195000, 315000} };			// Высоты выведения,м

	double DegreeToRadian(double degree) { return degree * M_PI / 180; }
	double RadianToDegree(double radian) { return radian * 180 / M_PI; }
// ПРАВЫЕ ЧАСТИ 
	double dx_dt(double Vx) { return Vx; }
	double dy_dt(double Vy) { return Vy; }
	double dVx_dt(double m, double t, double X, double Y) { return m_P / m * cos(tetta_control(t)) + get_gx(X, Y); }
	double dVy_dt(double m, double t, double X, double Y) { return m_P / m * sin(tetta_control(t)) + get_gy(X, Y); }
	double dm_dt(double betta) { return -betta; }
// Закон управления
	double tetta_control(double t);
// Вспомогательные функции
	double get_gx(double X, double Y) { return (-mu_M * X) / (pow(pow(X, 2) + pow(R_M + Y, 2), 1.5)); }
	double get_gy(double X, double Y) { return (-mu_M * (R_M + Y)) / (pow(pow(X, 2) + pow(R_M + Y, 2), 1.5)); }
	double get_V(double Vx, double Vy) { return sqrt(Vx * Vx + Vy * Vy); }
	double get_R(double X, double Y) { return sqrt(X * X + (R_M + Y) * (R_M + Y)); }
	double get_h(double X, double Y) { return sqrt(X * X + (R_M + Y) * (R_M + Y)) - R_M; }
	double get_Tetta(double params[]);
	double get_Tetta_C(double params[]) { return params[0] > 1e-8 ? atan(params[1] / params[0]) : M_PI_2; }
	double get_alpha(double t, double params[]) { return tetta_control(t) - get_Tetta_C(params); }
	double get_fi(double params[]) { return atan(params[2] / (params[3] + R_M)); }
	double get_nxa(double t, double params[]) { return m_P * cos(get_alpha(t, params)) / (params[4] * 1.623); }
	double get_nya(double t, double params[]) { return m_P * sin(get_alpha(t, params)) / (params[4] * 1.623); }
	double get_deltaV(double V, double V_orbit) { return V - V_orbit; }
	double get_deltaR(double R, double R_circle) { return R - R_circle; }
	double get_deltaTetta(double Tetta) { return Tetta - 0; }
	double get_deltaM(double m) { return m - m_m0 + m_mT_max; }	// Сколько осталось топлива,кг
// Интегрирование
	void Step(double& t, double dt, double params[], double old_params[]);	// Один шаг метода РК4
	struct Nevyazki		// Метод интегрирования имеет тип возвращаемого параметра - структура с невязками в конце траектории
	{
		double V;
		double r;
		double Tetta;
		double Massa;
	};
	Nevyazki RungeKutt(Hkon h_index, Pkon P_index);		// Интегрирование всей траектории полёта с одним набором параметров
	void StepFragment(double current, double old, double eps, bool isDiv, double params[], double old_params[],
		double &t, double &dt, double t_prev[], int &count, bool CutOff, int P_index);	// Дробление шага при различных условиях
// Обращение матрицы
	void InverseMatrix(double A[2][2], double A_1[2][2]);
// Решение Краевой задачи
	void SolveKZ();
// Решение задачи параметрической оптимизации алгоритмом Левенберга-Марквардта (МНК)
	void Optimization();
// Вычисление частных производных (в т.ч. вторых)
	double dg_dt1(double dt);
	double dg_dtPP(double dt);
	double d2g_dt12(double dt);
	double d2g_dtPP2(double dt);
	double d2g_dt1tPP(double dt);
	double d2g_dtPPt1(double dt);
// Запись в файлы
	void PrintReport(bool isFirstEnter, bool isExtr, bool isLast, double t, double params[]);		// Таблица в отчёт
	void PrintAllParams(bool isFirstEnter, bool isExtr, bool isLast, double t, double params[]);		// Все параметры
public:
	void Calculate();		// Основной метод класса (единственный public) - решение краевой задачи и оптимизация
};

