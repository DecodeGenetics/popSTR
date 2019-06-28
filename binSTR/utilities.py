

def tokenizer(s, c):
  i = 0
  while True:
    try:
      j = s.index(c, i)
    except ValueError:
      yield s[i:]
      return
    yield s[i:j]
    i = j + 1


def to_index(x, y):
  assert x <= y
  return x + (y + 1) * y // 2;


def bin2str(bin):
  return "+".join([str(x) for x in sorted(list(bin))])


def get_bin(bins, a):
  for i, bin in enumerate(bins):
    if a in bin:
      return i

  assert False
  return -1


def get_index2new_index(bins, n):
  index2new_index = {}

  for x in range(n):
    index2new_index[x] = get_bin(bins, x)

  return index2new_index


def get_gt(pl, n):
  i = 0
  for y in range(n):
    for x in range(y+1):
      if pl[i] == 0:
        return (x,y)

      i += 1

  assert False
  return (-1,-1)
