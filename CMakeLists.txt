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

file(GLOB project_sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB project_headers ${PROJECT_SOURCE_DIR}/include/*.hh)

file(GLOB alglib_sources /home/gabriel/Downloads/alglib-4.04.0.cpp.gpl/alglib-cpp/src/*.cc)
file(GLOB alglib_headers /home/gabriel/Downloads/alglib-4.04.0.cpp.gpl/alglib-cpp/src/*.h)


# Definindo o diretório do executável
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

# Criação do executável
add_executable(EttoreFit fit_comp.cpp ${project_sources} ${alglib_sources})

# Adicionando as flags manualmente
# Passando as opções de compilação corretamente para o CMake
string(REPLACE " " ";" ROOT_CFLAGS_LIST ${ROOT_CFLAGS})

# Incluindo os diretórios de header e as bibliotecas corretamente
target_include_directories(EttoreFit PRIVATE ${ROOT_CFLAGS_LIST} ${PROJECT_SOURCE_DIR}/include)
target_include_directories(EttoreFit PRIVATE /home/gabriel/Downloads/alglib-4.04.0.cpp.gpl/alglib-cpp/src)

target_compile_options(EttoreFit PRIVATE ${ROOT_CFLAGS_LIST})

# Linkando as bibliotecas do ROOT
target_link_libraries(EttoreFit PRIVATE ${ROOT_LIBS})




