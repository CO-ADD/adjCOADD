"""

"""
import logging
logger = logging.getLogger(__name__)


# ----------------------VLog 
class Validation_Log():

    def __init__(self,logProcess,logTypes= ['Error','Warning','Info']):
        self.logProcess = logProcess
        self.logTypes = logTypes
        self.nLogs = {}
        self.Logs  = {}
        self.Info  = {}

        for t in self.logTypes:
            self.nLogs[t] = 0
            self.Logs[t] = []
           

    def add_log(self, logType, logDesc, logItem, logHelp, ):
        lDict = {
            'Process': self.logProcess, 
            'Description': logDesc, 
            'Item': str(logItem), 
            'Help': logHelp,
#            'Time': datetime.now() 
            }
        logType = logType[0].upper()+logType[1:].lower()
        if logType in self.logTypes:
            self.Logs[logType].append(lDict)
            self.nLogs[logType] = self.nLogs[logType] + 1

    def select_unique(self,logTypes= ['Error','Warning', 'Info']):
        uLogs={}
        for t in logTypes:
            uLogs[t]=[]
            for l in self.Logs[t]:
                flAdd=False
                if len(uLogs[t])<1:
                    flAdd=True
                else:
                    flAdd=True
                    for u in uLogs[t]:
                        if u == l:
                            flAdd=False
                if flAdd:
                    uLogs[t].append(l)
        self.Logs = uLogs

    def reset(self):
        self.nLogs = {}
        self.Logs  = {}
        self.Info  = {}

        for t in self.logTypes:
            self.nLogs[t] = 0
            self.Logs[t] = []        

    def info(self,logTypes= ['Error','Warning', 'Info']):
        self.info={}
        for t in logTypes:
            self.info[t]=[]
            for l in self.Logs[t]:
                description=str(l['Description']).replace("'", "").replace('"', '')
                print_info=f"{l['Process']}_{description}_{l['Item']}_{l['Help']}"
                self.info[t].append(print_info) 

    def show(self,logTypes= ['Error','Warning', 'Info']):
        for t in logTypes:
            for l in self.Logs[t]:
                logger.info(f"{t:7s} - {l['Process']} : {l['Description']} ({l['Item']}) {l['Help']} ")
    
    def log_to_UI(self,logTypes= ['Error','Warning', 'Info']):
        info={} #info=[]
        for t in logTypes:
            # print(f"-- {t.upper():8} ({self.nLogs[t]:3}) ------------------------------------------------------")
            info[t]=[]
            for l in self.Logs[t]:
                print(f"{l['Process']}-{l['Description']} ({l['Item']}) {l['Help']} ")
                description=str(l['Description']).replace("'", "").replace('"', '')
                print_info=f"{l['Process']}_{description}_{l['Item']}_{l['Help']}"
                info[t].append(print_info) # info.append(print_info)
       
        return info
    