data = []
tmp = 0
with open('input.txt', 'rt') as inputfile:
        for inputline in inputfile:
                if tmp == 0:
                        tmp = 1
                else:
                        linedata = inputline.split()
                        data.append(float(linedata[2]))

print(data)
