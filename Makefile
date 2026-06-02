CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude
LDFLAGS = -lm

SRCS = src/sheaf.c src/agents.c src/consensus.c
OBJS = $(SRCS:.c=.o)

.PHONY: all test clean

all: test_sheaf_agents

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test_sheaf_agents: $(OBJS) tests/test_sheaf_agents.c
	$(CC) $(CFLAGS) tests/test_sheaf_agents.c $(OBJS) $(LDFLAGS) -o test_sheaf_agents

test: test_sheaf_agents
	./test_sheaf_agents

clean:
	rm -f $(OBJS) test_sheaf_agents
