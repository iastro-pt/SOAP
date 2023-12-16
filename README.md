# SOAP (Spot Oscillation And Planet) 2.0 Modernization Update

SOAP 2.0 has been modernized to enhance compatibility with contemporary software environments. The following README has been updated to reflect the new installation and execution processes that incorporate `pipenv` for Python dependency management, GSL (GNU Scientific Library) for advanced numerical analysis, and a transition to Python 3.x. Below you will find instructions that supersede the original guide wherever there is a conflict.

The original instructions by Xavier Dumusque, alongside historical and contextual references, can still be accessed for informational purposes at [the original SOAP repository](https://github.com/iastro-pt/SOAP/blob/ef4a54d41992654771aecf15b69af15337c999b3/README.md).

The goal of this fork is to eventually be worked into the upstream SOAP codebase.

---

## Introduction

SOAP 2.0 simulates the effects of stellar activity on radial velocity and photometric measurements, crucial for exoplanet detection and characterization. Refer to Dumusque et al. (2014, ApJ, 796, 132) for a detailed explanation of the scientific background and the code's operational principles.

## Prerequisites

Ensure you have a proper Git environment set up to clone the repository. You will also need Python 3.x and the GNU Scientific Library (GSL) installed on your system.

## Installation

### Download

Clone the current `master` branch of this repository or download a zip file:

```bash
git clone https://github.com/JonathanPorta/SOAP.git
```

Or download the zip file by navigating to `https://github.com/JonathanPorta/SOAP` and clicking the `Code` button followed by `Download ZIP`.

### Dependencies

SOAP 2.0 now uses `pipenv` for managing Python dependencies to ensure consistent environments for all users.

First, install `pipenv` if you haven't already:

```bash
pip install pipenv
```

Navigate to the repository directory and install the required Python dependencies:

```bash
cd SOAP
pipenv install
```

The required C libraries (GSL) should be installed on your system. For MacOS, use Homebrew:

```bash
brew install gsl
```

For Ubuntu and other Debian-based systems, use `apt-get`:

```bash
sudo apt-get install libgsl-dev
```

### Compilation

Activate the `pipenv` shell and compile the C code using the `setup.py` file:

```bash
pipenv shell
cd StarSpot
python setup.py build
cp build/lib.****/starspot.so .  # Replace **** with your operating system and architecture details
```

## Running SOAP 2.0

Configuration and usage instructions remain as detailed in the original README but ensure to use Python 3 syntax and functionalities.

## Outputs and Display

Outputs will be generated as described previously but consider running SOAP within the `pipenv` managed environment to avoid any library or version conflicts.

## Troubleshooting

Issues with GSL or Python libraries installation should be directed to system administrators or the respective library's support channels, as these are environment-specific and not issues with the SOAP codebase itself.

---

## Contact for the Modernized Version

For any issues or contributions related to the modernization of SOAP 2.0, please use the repository's Issues and Pull Requests sections on GitHub.

For the original SOAP code, contact Xavier Dumusque at the provided email address.
