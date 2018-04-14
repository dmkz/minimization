CXX ?= gcc
CXXFLAGS ?= -pthread -Ofast -std=c++14 -Wall -Wextra -fmax-errors=2 -I"include"

# Поддиректории, в которых также нужно осуществить сборку:
SUBDIRS := test-methods test-tms test-global

# Если в них нужно осуществлять сборку, то нужно и выполнять make clean
# Преобразуем цели, добавив суффикс, отвечающий за очистку:
CLEAN_SUFFIX := ___clean___
SUBDIRS_CLEAN := $(foreach dir,$(SUBDIRS),$(dir)$(CLEAN_SUFFIX))
EMPTY :=

.PHONY: all $(SUBDIRS) clean

all: $(SUBDIRS)
$(SUBDIRS):
		$(MAKE) -C $@
 
.c.o:
		$(CXX) $(CXXFLAGS) -c $< -o $@

clean: $(SUBDIRS_CLEAN)
$(SUBDIRS_CLEAN):
		$(MAKE) -C $(@:$(CLEAN_SUFFIX)=$(EMPTY)) clean