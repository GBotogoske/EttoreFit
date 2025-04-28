cmake_minimum_required(VERSION 3.10)
project(EttoreFit)

set(CMAKE_CXX_STANDARD 17)

# ROOT CFLAGS e LIBS
execute_process(
  COMMAND root-config --cflags
  OUTPUT_VARIABLE ROOT_CFLAGS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND root-config --libs
  OUTPUT_VARIABLE ROOT_LIBS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Definindo o diretório do executável
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

# Criação do executável
add_executable(EttoreFit fit_comp.cpp ${sources})

# Adicionando as flags manualmente
# Passando as opções de compilação corretamente para o CMake
string(REPLACE " " ";" ROOT_CFLAGS_LIST ${ROOT_CFLAGS})

# Incluindo os diretórios de header e as bibliotecas corretamente
target_include_directories(EttoreFit PRIVATE ${ROOT_CFLAGS_LIST} ${PROJECT_SOURCE_DIR}/include)
target_compile_options(EttoreFit PRIVATE ${ROOT_CFLAGS_LIST})

# Linkando as bibliotecas do ROOT
target_link_libraries(EttoreFit PRIVATE ${ROOT_LIBS})




