all: individualConstructions otherBins

CC = gcc
CFLAGS = -Wall -Wextra -ggdb -O2 -flto
LDFLAGS = -L/usr/local/lib -lflint -lm -Wl,-rpath,/usr/local/lib

OBJ_DIR = bin/objects
CONS_OBJ_DIR = $(OBJ_DIR)/cons
OBJ_SRC_DIR = src/core
CONS_OBJ_SRC_DIR = $(OBJ_SRC_DIR)/CFF_Constructions
INDIVIDUAL_CONSTRUCTION_BINS = bin/RS bin/STS bin/Sperner bin/ID bin/shortRS bin/Porat
OTHER_BINS = bin/getByN bin/getByT bin/writeTables bin/flint_test

MAIN_OBJS = $(OBJ_DIR)/CFF_Tables.o \
			$(OBJ_DIR)/CFF.o \
			$(OBJ_DIR)/CombinatorialAlgorithms.o \
			$(OBJ_DIR)/GeneralAlgorithms.o

CONS_OBJS = $(CONS_OBJ_DIR)/Additive_Construction.o \
			$(CONS_OBJ_DIR)/Doubling_Construction.o \
			$(CONS_OBJ_DIR)/Ext_By_One_Construction.o \
			$(CONS_OBJ_DIR)/Finite_Fields_Wrapper.o \
			$(CONS_OBJ_DIR)/Fixed_CFFs.o \
			$(CONS_OBJ_DIR)/ID_Construction.o \
			$(CONS_OBJ_DIR)/Kronecker_Construction.o \
			$(CONS_OBJ_DIR)/Optimized_Kronecker_Construction.o \
			$(CONS_OBJ_DIR)/Porat_Construction.o \
			$(CONS_OBJ_DIR)/Reed_Solomon_Construction.o \
			$(CONS_OBJ_DIR)/Short_Reed_Solomon_Construction.o \
			$(CONS_OBJ_DIR)/Sperner_Construction.o \
			$(CONS_OBJ_DIR)/STS_Construction.o

ALL_OBJS = $(MAIN_OBJS) $(CONS_OBJS)

$(OBJ_DIR)/%.o: $(OBJ_SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(CONS_OBJ_DIR)/%.o: $(CONS_OBJ_SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(INDIVIDUAL_CONSTRUCTION_BINS): bin/%: src/individualConstructions/%.c $(ALL_OBJS)
	$(CC) $(CFLAGS) $< $(ALL_OBJS) -o $@ $(LDFLAGS)

$(OTHER_BINS): bin/%: src/%.c $(ALL_OBJS)
	$(CC) $(CFLAGS) $< $(ALL_OBJS) -o $@ $(LDFLAGS)

constructionObjs: $(ALL_OBJS)
individualConstructions: $(INDIVIDUAL_CONSTRUCTION_BINS)
otherBins: $(OTHER_BINS)

clean:
	rm -f $(MAIN_OBJS) $(CONS_OBJS) $(INDIVIDUAL_CONSTRUCTION_BINS) $(OTHER_BINS)

clean_tables:
	rm -f tables/*

.PHONY: all clean clean_tables constructionObjs individualConstructions otherBins