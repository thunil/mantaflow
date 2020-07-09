#******************************************************************************
#
# MantaGen
# Copyright 2018 Steffen Wiewel, Moritz Becher, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
#******************************************************************************

# https://stackoverflow.com/questions/31875/is-there-a-simple-elegant-way-to-define-singletons/12850496#12850496
def singleton(cls):
    obj = cls()
    cls.__new__ = staticmethod(lambda cls: obj)
    try:
        del cls.__init__
    except AttributeError:
        pass
    return cls


from enum import Enum
class LogType(Enum):
    Info = 0
    Warning = 1
    Error = 2
    Fatal = 3

@singleton
class Logger(object):
    def __init__(self, mode=LogType.Info, summary=True):
        print("Initialized Logger")
        self.__mode = mode.value
        self.__summary = summary
        self.__messages = {}
        for entry in LogType:
            self.__messages[entry] = []
    def __del__(self):
        # print all warnings and errors on application end (depending on mode)
        if self.__summary:
            print("\nWarning/Error Summary:")
            for log_type in LogType:
                if log_type.value > self.__mode:
                    self.print_all(log_type)
    def __handle_log(self, log_type, msg):
        self.__messages[log_type].append(msg)
        if log_type.value >= self.__mode:
            format_msg = "[{}] {}".format(log_type.name, msg)
            if log_type == LogType.Fatal:
                assert False, format_msg
            else:
                print(format_msg)
    def info(self, msg):
        self.__handle_log(LogType.Info, msg)
    def warning(self, msg):
        self.__handle_log(LogType.Warning, msg)
    def error(self, msg):
        self.__handle_log(LogType.Error, msg)
    def fatal(self, msg):
        self.__handle_log(LogType.Fatal, msg)
    def print_all(self, log_type):
        messages = set(self.__messages[log_type])
        counts = [self.__messages[log_type].count(ele) for ele in messages]
        for count, msg in zip(counts, messages):
            print("[{}] {} [Count: {}]".format(log_type.name, msg, count))

def info(msg):
    Logger().info(msg)
def warning(msg):
    Logger().warning(msg)
def error(msg):
    Logger().error(msg)
def fatal(msg):
    Logger().fatal(msg)