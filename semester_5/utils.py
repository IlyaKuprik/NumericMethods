def read_data(path):
    data = [[], []]
    f = open(path, 'r')
    print(path)
    for i in range(1000):
        line = f.readline().strip().split(" ")
        if line[0] == '' or line[1] is None:
            continue
        data[0].append(float(line[0]))
        data[1].append(float(line[1]))
    f.close()
    return data
