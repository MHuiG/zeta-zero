
f = open("./results/test.md", 'w')
f.write("|  No.   | Zero  |\n|  ----  | ----  |\n")

def write_zero(no,zero):
  f.write("|  "+str(no)+" | 1/2+"+str(zero)+"i |\n")


print("hello")
write_zero(1,11)
write_zero(2,22)
