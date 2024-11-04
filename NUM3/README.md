# NUM3

Program `NUM3` służy do rozwiązywania zadania z algebry liniowej, w którym wyznaczamy wektor `y = A^{-1} x` dla trójdiagonalnej macierzy `A` oraz obliczamy wyznacznik macierzy `A`.

## Wymagania

- Kompilator `g++` wspierający standard C++11.
- System operacyjny z obsługą `make` (Linux, macOS, lub WSL na Windowsie).

## Struktura katalogów

- `NUM3.cpp` - plik źródłowy zawierający kod programu.
- `build/` - folder, w którym zostanie umieszczony skompilowany plik wykonywalny.
- `Makefile` - plik make umożliwiający kompilację i uruchomienie programu.
  
## Instrukcje

### 1. **Kompilacja programu:**

Aby skompilować program, w terminalu wpisz:

```bash
make
```

Spowoduje to skompilowanie pliku `NUM3.cpp` i umieszczenie pliku wykonywalnego `NUM3` w folderze `build`.

### 2. **Kompilacja i uruchomienie programu:**

Aby skompilować i natychmiast uruchomić program, użyj:

```bash
make run
```

Program zostanie skompilowany (jeśli to konieczne), a następnie uruchomiony. Wynik zostanie wyświetlony w terminalu.

### 3. **Czyszczenie plików wynikowych:**

Aby usunąć wszystkie wygenerowane pliki i folder `build`, wpisz:

```bash
make clean
```

Spowoduje to usunięcie folderu `build` wraz ze skompilowanym plikiem wykonywalnym.
