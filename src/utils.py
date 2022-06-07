
def each_cons(xs, n):
    return [xs[i:i+n] for i in range(len(xs)-n+1)]