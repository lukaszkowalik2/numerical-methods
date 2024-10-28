# Generowanie Wektora Zaburzenia ∆b

Ten projekt w języku C++ generuje losowy wektor zaburzenia `∆b` o małej normie euklidesowej (ok. \(10^{-6}\)) przy użyciu biblioteki **Eigen**.

## Struktura Projektu

project_root/
├── lib/
│   └── eigen-3.4.0/     # Folder zawierający bibliotekę Eigen
├── build/               # Folder z plikami wynikowymi (tworzony automatycznie)
├── NUM2.cpp             # Plik źródłowy z kodem programu
├── Makefile             # Plik Makefile do kompilacji i uruchomienia programu
└── README.md            # Dokumentacja projektu

## Wymagania

- **Kompilator C++** z obsługą standardu C++11 lub nowszego (np. `g++`).
- **Biblioteka Eigen** (dostarczona w projekcie w katalogu `lib/eigen-3.4.0`).

## Instrukcje Kompilacji i Uruchomienia

W katalogu głównym projektu możesz użyć poniższych poleceń do kompilacji, uruchomienia i czyszczenia projektu:

### 1. Kompilacja Programu

Aby skompilować program, użyj:

```bash
make
```

Polecenie to utworzy katalog `build` i wygeneruje plik wykonywalny `NUM2` w folderze `build`.

### 2. Uruchomienie Programu

Aby skompilować i natychmiast uruchomić program, użyj:

```bash
make run
```

### 3. Czyszczenie Plików Wynikowych

Aby usunąć skompilowane pliki z katalogu `build`, użyj:

```bash
make clean
```
