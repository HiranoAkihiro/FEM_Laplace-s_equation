FC = gfortran
FFLAGS = -O2 -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow

MOD_DIR = -J ./include
BIN_DIR = ./bin
SRC_DIR = ./src
OBJ_DIR = ./obj
BIN_LIST = poisson_eq
TARGET = $(addprefix $(BIN_DIR)/, $(BIN_LIST))

SRC_LIST = \
utils.f90 \
io.f90 \
matrix.f90 \
valid.f90 \
main_poisson.f90 \

SOURCES = $(addprefix $(SRC_DIR)/, $(SRC_LIST))
OBJS = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES:.f90=.o))

$(TARGET): $(OBJS)
	$(FC) -o $@ $(OBJS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(MOD_DIR) -o $@ -c $<

clean:
	rm $(OBJS) $(TARGET) ./include/*.mod

.PHONY: clean