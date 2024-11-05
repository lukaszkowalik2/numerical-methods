# NUM3

Program `NUM3` służy do rozwiązywania zadania z algebry liniowej, w którym wyznaczamy wektor \( y = A^{-1} x \) dla macierzy czterodiagonalnej \( A \) oraz obliczamy wyznacznik tej macierzy. Program umożliwia porównanie różnych algorytmów numerycznych i mierzy ich czas działania.

## Wymagania

### Python

- `matplotlib` do generowania wykresów
- `pandas` do analizy danych czasów działania

Aby zainstalować wymagane pakiety Python, użyj:

```bash
pip install matplotlib pandas
```

## Struktura katalogów

- `NUM3.cpp` - plik źródłowy obliczający wektor \( y = A^{-1} x \) i wyznacznik dla \( N = 300 \).
- `NUM3_timed.cpp` - plik źródłowy, który wykonuje obliczenia dla różnych wartości \( N \), generując plik `times.csv`.
- `algorithms.h` - plik nagłówkowy zawierający implementacje wszystkich algorytmów (eliminacja Gaussa, faktoryzacja LU, zmodyfikowany algorytm Thomasa).
- `lib/eigen-3.4.0/` - folder z biblioteką Eigen, potrzebną do obsługi algebry liniowej.
- `build/` - folder, w którym zostanie umieszczony skompilowany plik wykonywalny.
- `Makefile` - plik make umożliwiający kompilację i uruchomienie programu.
- `generate_plot.py` - skrypt Python generujący wykres zależności czasowych dla różnych wartości \( N \), korzystający z danych zapisanych w pliku `times.csv`.

## Instrukcje

### 1. Kompilacja programu

Aby skompilować program, w terminalu wpisz:

```bash
make
```

Spowoduje to skompilowanie plików `NUM3.cpp` i `NUM3_timed.cpp` oraz umieszczenie plików wykonywalnych w folderze `build`.

### 2. Kompilacja i uruchomienie programu dla \( N = 300 \)

Aby skompilować i natychmiast uruchomić program dla \( N = 300 \), użyj:

```bash
make run
```

Program `NUM3.cpp` zostanie skompilowany (jeśli to konieczne) i uruchomiony. Wyniki dla \( N = 300 \) zostaną wyświetlone w terminalu.

### 3. Generowanie pliku `times.csv` i wykresu

Aby wykonać obliczenia dla różnych wartości \( N \) (co 10 od \( N = 10 \) do \( N = 300 \)) i wygenerować wykres, użyj:

```bash
make chart
```

To polecenie uruchomi program `NUM3_timed.cpp`, który zapisze czasy działania dla poszczególnych algorytmów w pliku `times.csv`. Następnie skrypt `generate_plot.py` wygeneruje wykres zależności czasowych dla różnych algorytmów.

### 4. Czyszczenie plików wynikowych

Aby usunąć wszystkie wygenerowane pliki i folder `build`, wpisz:

```bash
make clean
```

Spowoduje to usunięcie folderu `build` wraz ze skompilowanymi plikami wykonywalnymi oraz plikami wynikowymi.
