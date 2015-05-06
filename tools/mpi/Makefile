## Preamble {{{
## =============================================================================
##       @file Makefile
## -----------------------------------------------------------------------------
##     @author Alexander Shukaev <"alexander" "." "shukaev" "@" "tum" "." "de">
## -----------------------------------------------------------------------------
## @maintainer Alexander Shukaev <"alexander" "." "shukaev" "@" "tum" "." "de">
## -----------------------------------------------------------------------------
##  @copyright Copyright (C) 2015,
##             Alexander Shukaev <"alexander" "." "shukaev" "@" "tum" "." "de">.
##             All rights reserved.
## -----------------------------------------------------------------------------
##    @license This program is free software: you can redistribute it and/or
##             modify it under the terms of the GNU General Public License as
##             published by the Free Software Foundation, either version 3 of
##             the License, or (at your option) any later version.
##
##             This program is distributed in the hope that it will be useful,
##             but WITHOUT ANY WARRANTY; without even the implied warranty of
##             MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
##             General Public License for more details.
##
##             You should have received a copy of the GNU General Public License
##             along with this program.  If not, see
##             <http://www.gnu.org/licenses/>.
## =============================================================================
## }}} Preamble

# Directories {{{
# ==============================================================================
BUILD_DIR   := .build
INSTALL_DIR := .install
SOURCE_DIR  := .
# ==============================================================================
# }}} Directories

# Compilers {{{
# ==============================================================================
ifndef CXX
  CXX := $(shell which mpicxx 2> /dev/null)
endif

ifndef CXX
  CXX := $(shell which mpicc 2> /dev/null)
endif

ifndef CXX
  CXX := $(shell which icpc 2> /dev/null)
endif

ifndef CXX
  CXX := $(shell which icc 2> /dev/null)
endif

ifndef CXX
  CXX := $(shell which g++ 2> /dev/null)
endif

ifndef CXX
  CXX := $(shell which gcc 2> /dev/null)
endif

ifndef CXX
  $(error "C++ compiler not found")
endif
# ==============================================================================
# }}} Compilers

# Flags {{{
# ==============================================================================
ifndef CXXFLAGS
  CXXFLAGS :=                                                                  \
    -march=native                                                              \
    -O3                                                                        \
    -fno-rtti                                                                  \
    -fomit-frame-pointer                                                       \
    -funroll-loops
endif

override                                                                       \
CXXFLAGS :=                                                                    \
  $(CXXFLAGS)                                                                  \
  -std=c++11                                                                   \
  -pipe                                                                        \
  -Wall                                                                        \
  -Wextra                                                                      \
  -Wpedantic                                                                   \
  -Werror

CPPFLAGS := $(CPPFLAGS) -DMPICH_SKIP_MPICXX

LDFLAGS := $(LDFLAGS) $(CXXFLAGS)
# ==============================================================================
# }}} Flags

# Configurations {{{
# ==============================================================================
# `acceptor' {{{
# ==============================================================================
ACCEPTOR := acceptor

ACCEPTOR_SOURCES := main.cpp

ACCEPTOR_BUILD_DIR  := $(BUILD_DIR)/$(ACCEPTOR)
ACCEPTOR_SOURCE_DIR := $(SOURCE_DIR)/$(ACCEPTOR)

ACCEPTOR_OBJECTS := $(addsuffix .o,                      $(ACCEPTOR_SOURCES))
ACCEPTOR_OBJECTS := $(addprefix $(ACCEPTOR_BUILD_DIR)/,  $(ACCEPTOR_OBJECTS))
ACCEPTOR_SOURCES := $(addprefix $(ACCEPTOR_SOURCE_DIR)/, $(ACCEPTOR_SOURCES))
ACCEPTOR         := $(addprefix $(ACCEPTOR_BUILD_DIR)/,  $(ACCEPTOR))

ACCEPTOR_CXXFLAGS := $(CXXFLAGS)
ACCEPTOR_CPPFLAGS := $(CPPFLAGS)
ACCEPTOR_LDFLAGS  := $(LDFLAGS)
# ==============================================================================
# }}} `acceptor'

# `requester' {{{
# ==============================================================================
REQUESTER := requester

REQUESTER_SOURCES := main.cpp

REQUESTER_BUILD_DIR  := $(BUILD_DIR)/$(REQUESTER)
REQUESTER_SOURCE_DIR := $(SOURCE_DIR)/$(REQUESTER)

REQUESTER_OBJECTS := $(addsuffix .o,                       $(REQUESTER_SOURCES))
REQUESTER_OBJECTS := $(addprefix $(REQUESTER_BUILD_DIR)/,  $(REQUESTER_OBJECTS))
REQUESTER_SOURCES := $(addprefix $(REQUESTER_SOURCE_DIR)/, $(REQUESTER_SOURCES))
REQUESTER         := $(addprefix $(REQUESTER_BUILD_DIR)/,  $(REQUESTER))

REQUESTER_CXXFLAGS := $(CXXFLAGS)
REQUESTER_CPPFLAGS := $(CPPFLAGS)
REQUESTER_LDFLAGS  := $(LDFLAGS)
# ==============================================================================
# }}} `requester'
# ==============================================================================
# }}} Configurations

# Targets {{{
# ==============================================================================
.PHONY: all
all: acceptor requester

.PHONY: clean
clean:
	@rm -f -r -d $(BUILD_DIR)/

.PHONY: install
install: $(INSTALL_DIR) all
	@install -m 0755 $(ACCEPTOR) $(REQUESTER) $(INSTALL_DIR)/

.PHONY: uninstall
uninstall:
	@rm -f -r -d $(INSTALL_DIR)/

.PHONY: acceptor
acceptor: $(ACCEPTOR)

.PHONY: requester
requester: $(REQUESTER)

$(ACCEPTOR_BUILD_DIR):
	@mkdir -p $(ACCEPTOR_BUILD_DIR)/

$(REQUESTER_BUILD_DIR):
	@mkdir -p $(REQUESTER_BUILD_DIR)/

$(INSTALL_DIR):
	@mkdir -p $(INSTALL_DIR)/

$(ACCEPTOR_BUILD_DIR)/%.cpp.o: $(ACCEPTOR_SOURCE_DIR)/%.cpp
	$(CXX) $(ACCEPTOR_CXXFLAGS) $(ACCEPTOR_CPPFLAGS) -c $< -o $@

$(REQUESTER_BUILD_DIR)/%.cpp.o: $(REQUESTER_SOURCE_DIR)/%.cpp
	$(CXX) $(REQUESTER_CXXFLAGS) $(REQUESTER_CPPFLAGS) -c $< -o $@

$(ACCEPTOR): $(ACCEPTOR_BUILD_DIR) $(ACCEPTOR_OBJECTS)
	$(CXX) -o $@ $(ACCEPTOR_OBJECTS) $(ACCEPTOR_LDFLAGS)

$(REQUESTER): $(REQUESTER_BUILD_DIR) $(REQUESTER_OBJECTS)
	$(CXX) -o $@ $(REQUESTER_OBJECTS) $(REQUESTER_LDFLAGS)
# ==============================================================================
# }}} Targets
