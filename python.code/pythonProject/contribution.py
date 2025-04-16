#参数a是一个元组，如a = (6,5,4,3,2,1)
def con(a):
    b = [0]  len(a)
    b [ 0 ] = a [ 0 ] 2sum(a)2
    for i in range(len(a))
        if i!= 0
            m = a[ i+1]
            mm = a[ i ]
            b[i] = (sum(m)  2 - sum(mm)  2)  sum(a)  2
    print(b)
    return b