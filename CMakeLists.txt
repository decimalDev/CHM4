# CMakeList.txt: проект CMake для lab4; включите исходный код и определения,
# укажите здесь логику для конкретного проекта.
#
cmake_minimum_required (VERSION 3.8)

project ("lab4")

# Добавьте источник в исполняемый файл этого проекта.
add_executable (lab4 "lab4.cpp" "lab4.h" "Task I.cpp" "Task II.cpp" "Task III.cpp" "Matrix algebra.cpp")

# TODO: Добавьте тесты и целевые объекты, если это необходимо.
