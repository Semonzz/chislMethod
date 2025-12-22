#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

using namespace std;

void sysout(ofstream& fout, const vector<vector<double>>& a, const vector<double>& y)
{
    int n = a.size();
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fout << a[i][j] << "*x" << j;
            if (j < n - 1)
                fout << " + ";
        }
        fout << " = " << y[i] << '\n';
    }
    return;
}

void sysout_console(const vector<vector<double>>& a, const vector<double>& y)
{
    int n = a.size();
    for (int i = 0; i < n; i++)
    {
        cout << "Уравнение " << (i + 1) << ": ";
        for (int j = 0; j < n; j++)
        {
            cout << setw(6) << a[i][j] << " * x" << (j + 1);
            if (j < n - 1)
                cout << " + ";
        }
        cout << " = " << setw(6) << y[i] << "\n";
    }
    cout << "\n";
    return;
}

void printMatrix(const vector<vector<double>>& mat, const vector<double>& rhs) {
    int n = mat.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            cout << setw(8) << fixed << setprecision(4) << mat[i][j];
        cout << " | " << setw(8) << rhs[i] << "\n";
    }
    cout << "------------------------\n";
}

vector<double> gauss(vector<vector<double>> matrix, vector<double> rhs, bool showSteps = true) {
    int n = matrix.size();
    const double eps = 1e-12;

    if (showSteps) {
        cout << "\n[Гаусс] Исходная система:\n";
        printMatrix(matrix, rhs);
    }

    // Прямой ход: превращаем матрицу в верхнетреугольную
    for (int col = 0; col < n; col++) {
        // --- Шаг 1: Найти строку с максимальным элементом в текущем столбце (частичный выбор главного элемента)
        int bestRow = col;
        double maxInColumn = abs(matrix[col][col]);
        for (int row = col + 1; row < n; row++) {
            if (abs(matrix[row][col]) > maxInColumn) {
                maxInColumn = abs(matrix[row][col]);
                bestRow = row;
            }
        }

        // --- Шаг 2: Проверка вырожденности
        if (maxInColumn < eps) {
            bool hasNonZeroRHS = false;
            for (int row = col; row < n; row++) {
                if (abs(rhs[row]) > eps) {
                    hasNonZeroRHS = true;
                    break;
                }
            }
            if (hasNonZeroRHS)
                cout << "Система несовместна! Решений нет.\n";
            else
                cout << "Система имеет бесконечно много решений.\n";
            return {};
        }

        // --- Шаг 3: Переставить строки, если нужно
        if (bestRow != col) {
            swap(matrix[col], matrix[bestRow]);
            swap(rhs[col], rhs[bestRow]);
            if (showSteps) {
                cout << "[Гаусс] Перестановка строк " << col + 1 << " и " << bestRow + 1 << "\n";
                printMatrix(matrix, rhs);
            }
        }

        // --- Шаг 4: Сделать ведущий элемент = 1 (нормализация)
        double pivot = matrix[col][col];
        for (int j = col; j < n; j++) {
            matrix[col][j] /= pivot;
        }
        rhs[col] /= pivot;

        if (showSteps) {
            cout << "[Гаусс] Нормализация строки " << col + 1 << " (pivot = " << pivot << ")\n";
            printMatrix(matrix, rhs);
        }

        // --- Шаг 5: Обнулить текущий столбец во всех строках ниже
        for (int row = col + 1; row < n; row++) {
            double multiplier = matrix[row][col]; // сколько раз строка col "влезает" в строку row
            if (abs(multiplier) < eps) continue;  // уже ноль — пропускаем

            for (int j = col; j < n; j++) {
                matrix[row][j] -= multiplier * matrix[col][j];
            }
            rhs[row] -= multiplier * rhs[col];

            if (showSteps) {
                cout << "[Гаусс] Обнулили столбец " << col + 1
                    << " в строке " << row + 1
                    << " (множитель = " << multiplier << ")\n";
                printMatrix(matrix, rhs);
            }
        }
    }

    // Обратный ход: находим решение снизу вверх
    vector<double> solution(n);
    for (int row = n - 1; row >= 0; row--) {
        solution[row] = rhs[row];
        for (int col = row + 1; col < n; col++) {
            solution[row] -= matrix[row][col] * solution[col];
        }
        if (showSteps) {
            cout << "x" << row + 1 << " = " << fixed << setprecision(6) << solution[row] << "\n";
        }
    }

    return solution;
}

vector<double> jacobi(const vector<vector<double>>& mat, const vector<double>& rhs,
    double tol = 1e-10, int maxIter = 1000, bool show = true) {
    int n = mat.size();
    const double eps = 1e-12;

    // Проверка: диагональ не должна быть нулём
    for (int i = 0; i < n; i++) {
        if (abs(mat[i][i]) < eps) {
            cout << "Якоби невозможен: нулевая диагональ в строке " << i << "\n";
            return {};
        }
    }

    vector<double> x(n, 0.0);  // начальное приближение

    if (show) {
        cout << "\n[Якоби] Старт: все x = 0\n";
    }

    for (int iter = 1; iter <= maxIter; iter++) {
        vector<double> x_next(n);
        double max_err = 0.0;

        // Считаем новое приближение
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) sum += mat[i][j] * x[j];
            }
            x_next[i] = (rhs[i] - sum) / mat[i][i];
            max_err = max(max_err, abs(x_next[i] - x[i]));
        }

        // Проверка сходимости
        if (max_err < tol) {
            if (show) cout << "[Якоби] Сошёлся на итерации " << iter << "\n";
            return x_next;
        }

        // Вывод каждые 5 итераций
        if (show && (iter <= 5 || iter % 5 == 0)) {
            cout << "Итерация " << iter << ": ";
            for (int i = 0; i < n; i++)
                cout << "x" << i + 1 << "=" << fixed << setprecision(6) << x_next[i] << " ";
            cout << "| ошибка = " << scientific << setprecision(2) << max_err << "\n";
        }

        x = x_next;
    }

    if (show) cout << "[Якоби] Достигнут лимит итераций (" << maxIter << ")\n";
    return x;
}

int chooseMethod() {
    int c;
    cout << "\n1. Гаусс\n2. Якоби\n3. Сравнить\nВыбор (1-3): ";
    while (!(cin >> c) || c < 1 || c > 3) {
        cout << "Неверно. Введите 1, 2 или 3: ";
        cin.clear();
        cin.ignore(1000, '\n');
    }
    return c;
}

int main() {
    system("chcp 1251 > nul");

    ifstream fin("input.txt");
    if (!fin) {
        cout << "Нет input.txt\n";
        return 1;
    }

    int n;
    fin >> n;
    if (n <= 0) {
        cout << "Некорректная размерность\n";
        return 1;
    }

    vector<vector<double>> a(n, vector<double>(n));
    vector<double> y(n);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            fin >> a[i][j];
    for (int i = 0; i < n; i++)
        fin >> y[i];
    fin.close();

    cout << "Система " << n << "x" << n << ":\n";
    sysout_console(a, y);

    int method = chooseMethod();
    vector<double> x;

    ofstream fout("output.txt");

    if (method == 1) {
        cout << "\n--- Гаусс ---\n";
        fout << "Метод: Гаусс\n";
        x = gauss(a, y, true);
    }
    else if (method == 2) {
        cout << "\n--- Якоби ---\n";
        fout << "Метод: Якоби\n";
        x = jacobi(a, y, 1e-10, 1000, true);
    }
    else {
        cout << "\n--- Сравнение ---\n";
        fout << "Сравнение методов\n";

        vector<double> x1 = gauss(a, y, false);
        vector<double> x2 = jacobi(a, y, 1e-10, 1000, false);

        if (!x1.empty() && !x2.empty()) {
            cout << "\nСравнение:\n";
            cout << "      Гаусс     Якоби\n";
            for (int i = 0; i < n; i++) {
                cout << "x" << i + 1 << ": "
                    << fixed << setprecision(6)
                    << setw(10) << x1[i]
                    << setw(10) << x2[i] << "\n";
                fout << "x" << i + 1 << ": Гаусс=" << x1[i] << ", Якоби=" << x2[i] << "\n";
            }
            x = x1;
        }
        else {
            cout << "Один из методов не сработал.\n";
            fout << "Ошибка в одном из методов.\n";
        }
    }

    if (!x.empty()) {
        cout << "\nРешение:\n";
        for (int i = 0; i < n; i++) {
            cout << "  x" << i + 1 << " = " << fixed << setprecision(6) << x[i] << "\n";
            fout << "x[" << i << "] = " << x[i] << "\n";
        }
    }
    else {
        cout << "Решение не найдено.\n";
        fout << "Решение не найдено.\n";
    }

    fout.close();
    cout << "\nРезультаты в output.txt\n";
    cin.get(); cin.get();
    return 0;
}