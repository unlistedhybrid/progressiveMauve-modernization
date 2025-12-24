# progressiveMauve (Modernized v2.5.0)

**A modernized, deterministic, and cross-platform version of the progressiveMauve multiple genome aligner.**

This repository hosts version 2.5.0 of **progressiveMauve**, originally developed by Aaron Darling at the [Darling Lab](https://darlinglab.org/mauve/mauve.html). This modernization project updates the legacy 2015 codebase to C++17 standards, ensuring it compiles on modern systems while guaranteeing bit-exact reproducibility across all platforms and CPU architectures.

## 🚀 Key Improvements in v2.5.0

This release focuses on stability, security, and deterministic output. It is fully compatible with existing Mauve workflows but robust enough for modern high-performance environments.

* **Bit-Exact Reproducibility:** Alignment results are now 100% deterministic. You will get the exact same alignment output regardless of whether you run the tool on Windows, macOS (Intel/Apple Silicon), or Linux.
* **Modern Compatibility:** The codebase has been updated to **C++17**, fixing build errors on modern GCC, Clang, and MSVC compilers.
* **Single-File Binary:** All dependencies (including Boost 1.89) are statically linked. There are no external libraries to install; just download the executable and run.
* **Crash Prevention:** Replaced legacy static buffers with dynamic resizing to handle arbitrarily long sequences without crashing. Fixed "cached file" errors by enforcing fresh data calculation on every run.
* **Security Fixes:** Replaced unsafe memory functions (like `sprintf`) with secure alternatives to prevent buffer overflows.

## 📥 Installation

Go to the [**Releases Page**](https://github.com/unlistedhybrid/progressiveMauve-modernization/releases) and download the archive for your operating system.

* **Windows:** Extract the `.zip` file. Run `progressiveMauve.exe` from the Command Prompt or PowerShell.
* **macOS / Linux:** Extract the `.tar.gz` file. The binary is pre-compiled and ready to run (no `chmod` required).

## 📖 Usage

Basic usage remains identical to the original version:

```bash
./progressiveMauve --output=my_alignment.xmfa genome1.fasta genome2.fasta

```

For a full list of options, run:

```bash
./progressiveMauve --help

```

## 📜 Attribution & License

**progressiveMauve** is free software released under the **GNU General Public License v2.0 (GPLv2)**.

This project is a modernization of the original work by **Aaron Darling** and colleagues.

* **Original Project Website:** [Darling Lab](https://darlinglab.org/mauve/mauve.html)
* **Original Source Code:** [SourceForge](https://sourceforge.net/projects/mauve/)
* **Original GitHub:** [koadman/mauveAligner](https://github.com/koadman/mauveAligner)

*The `muscle` directory contains a public domain library refactorization of Robert Edgar's MUSCLE v3.7 program.*

## 🤝 Contributing & Legacy

This repository aims to keep progressiveMauve functional for the future. If you are a developer looking to contribute, please see the source code, which has been refactored for modern C++ standards.

*This project is maintained as a modernized fork to ensure continued accessibility on newer hardware.*
