CC     = gcc
CFLAGS = -O0 -Iincludes
LFLAGS = -lm

PROG   = labZero
SRCS   = src/utils.c src/ZeroFuncao.c
OBJS   = $(SRCS:.c=.o)
MAIN   = labZero.o

DISTFILES = src/*.c includes/*.h Makefile assets/polinomios.dat
DISTDIR = labZero_dist

.PHONY: clean purge distclean all

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

$(PROG): $(OBJS) $(MAIN)
	$(CC) -o $@ $^ $(LFLAGS)

clean:
	@rm -f *~ *.bak src/*~ includes/*~

purge: clean
	@rm -f src/*.o *.o core a.out $(PROG)

dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tar.gz) ..."
	@mkdir -p $(DISTDIR)
	@cp -r $(DISTFILES) $(DISTDIR)/
	@tar -czvf $(DISTDIR).tar.gz $(DISTDIR)
	@rm -rf $(DISTDIR)
