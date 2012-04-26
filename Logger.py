import sys
class Logger(object):
    def __init__(self, filename="Default.log"):
        self.log = open(filename, "w")
        
    def write(self, message):
        sys.__stdout__.write(message)
        self.log.write(message)

    def close(self):
        self.log.close()
        sys.stdout = sys.__stdout__
