<!--
 !-- Copyright (C) 2019  Jimmy Aguilar Mena
 !--
 !-- This program is free software: you can redistribute it and/or modify
 !-- it under the terms of the GNU General Public License as published by
 !-- the Free Software Foundation, either version 3 of the License, or
 !-- (at your option) any later version.
 !--
 !-- This program is distributed in the hope that it will be useful,
 !-- but WITHOUT ANY WARRANTY; without even the implied warranty of
 !-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 !-- GNU General Public License for more details.
 !--
 !-- You should have received a copy of the GNU General Public License
 !-- along with this program.  If not, see <http://www.gnu.org/licenses/>.
  -->

# Readme


This is the basic centralized set of benchmarks for
[OmpSs-2@Clusters](https://github.com/bsc-pm/nanos6/blob/master/docs/cluster/README-CLUSTER.md).

## Installation

Right now it requires a set of dependencies and requirements to
work. We will try to reduce the dependencies as much as possible in
the future.

1. [OmpSs-2@Clusters
requirements](https://github.com/bsc-pm/nanos6/blob/master/docs/cluster/README-CLUSTER.md#system-requirements).

2. [Mercurium](https://pm.bsc.es/mcxx)

3. [Nanos6](https://pm.bsc.es/ftp/ompss-2/doc/user-guide/build/nanos6.html)

To build the code you need first to install all the dependencies above. Then just:

```console
git clone --recursive 'url_to_this_repo'
cd nanos6-cluster-benchmarks
mkdir build
cd build
cmake ..
make
```

When cloning `--recursive` some submodules will be downloaded and compiled automatically.

4. [ArgParserC](https://github.com/Ergus/ArgParserC.git)

5. [cmacros](https://github.com/Ergus/cmacros.git)

6. [ArgParserBash](https://github.com/Ergus/ArgParseBash.git)

After this all the benchmarks should be build inside a directory with
the same name than the original one in the project's root directory.

## Adding benchmarks

To add a new benchmark you only need to create a new directory inside
the project's root directory and add a CMakeLists.txt within it with
the build instructions.

You don't need to modify anything in the main CMakeLists.txt or
anywhere else outside that directory. If you want to add similar
benchmarks that reuses some code between them you can add them in the
same directory.

### Contribution hints

1. Please try to rely in the **utils** functionalities when possible
   instead of using external no standard tools.

2. Make the benchmarks keeping readability and simplicity as much as
   possible (read [The Zen of
   Python](https://www.python.org/dev/peps/pep-0020/) or execute
   `python -c 'import this'`).
   
3. You can use C (.h .c) and C++ (.hpp .cpp) mainly, please respect
   the filenames.
   
4. Use the [linux kernel coding
   style](https://www.kernel.org/doc/html/v4.10/process/coding-style.html).

## Tests with CTest

Some very simple set of tests are added when a **Python** and
**NumPy** are detected. These are very helpfull during developement
and the addition is optional.

It is recommended to add such tests in order to check for correctness
in the benchmarks, the runtime and the executions during
modifications.


