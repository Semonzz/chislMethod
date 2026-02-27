#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double f1(double x, double y) {
    return sin(x + y) - 1.1 * x - 0.1;
}

double f2(double x, double y) {
    return x * x + y * y - 1;
}

bool newton(double& x, double& y, int maxIter = 50, double eps = 1e-6) {
    cout << "Iter |     x      |     y      |    dx      |    dy\n";
    cout << "-----|------------|------------|------------|------------\n";

    for (int i = 0; i < maxIter; i++) {
        // Невязки
        double F1 = f1(x, y);
        double F2 = f2(x, y);

        // Матрица Якоби
        double J11 = cos(x + y) - 1.1;  // ∂f1/∂x
        double J12 = cos(x + y);         // ∂f1/∂y
        double J21 = 2 * x;                // ∂f2/∂x
        double J22 = 2 * y;                // ∂f2/∂y

        // Определитель
        double det = J11 * J22 - J12 * J21;
        if (fabs(det) < 1e-12) return false;

        // Обратная матрица
        double invJ11 = J22 / det;
        double invJ12 = -J12 / det;
        double invJ21 = -J21 / det;
        double invJ22 = J11 / det;

        // Приращения
        double dx = -(invJ11 * F1 + invJ12 * F2);
        double dy = -(invJ21 * F1 + invJ22 * F2);

        // Вывод итерации
        cout << setw(4) << i + 1 << " | "
            << setw(10) << fixed << setprecision(6) << x << " | "
            << setw(10) << y << " | "
            << setw(10) << setprecision(6) << dx << " | "
            << setw(10) << dy << "\n";

        // Обновление
        x += dx;
        y += dy;

        // Проверка сходимости
        if (sqrt(dx * dx + dy * dy) < eps) {
            cout << "-----|------------|------------|------------|------------\n";
            cout << "-> Сходимость за " << i + 1 << " итераций\n";
            return true;
        }

        if (fabs(dx) > 1e6 || fabs(dy) > 1e6 || isnan(x + y))
            return false;
    }
    return false;
}

int main() {
    setlocale(LC_ALL, "ru_RU");
    double x = 0.5, y = 0.5;

    cout << "Система: sin(x+y) - 1.1x = 0.1\n";
    cout << "         x^2 + y^2 = 1\n\n";
    cout << "Start: x=" << x << ", y=" << y << "\n\n";

    if (newton(x, y)) {
        cout << fixed << setprecision(8);
        cout << "\n=== РЕШЕНИЕ ===\n";
        cout << "x = " << x << "\n";
        cout << "y = " << y << "\n\n";

        cout << "Проверка:\n";
        cout << "f1 = " << f1(x, y) << " (должно быть 0)\n";
        cout << "f2 = " << f2(x, y) << " (должно быть 0)\n";
    }
    else {
        cout << "\nНЕ СОШЛОСЬ!\n";
    }

    return 0;
}