"""

"""
import pandas as pd
import logging
logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------- 
class Validation_Log():
    """
    In-process Logging class to capture outcomes of validation and processing tasks
        as logTypes = ['Error','Warning','Info']
    """
    
    LOG_FIELDS = ['Process','Text','Item','Note','Help']
    LOG_TYPES = ['Error','Warning','Info']
# ---------------------------------------------------------------------------

    #-----------------------------------------------------
    # Inits the log with a logProcess as Name
    #-----------------------------------------------------
    def __init__(self,logProcess,logTypes= LOG_TYPES):
        self.logProcess = logProcess
        self.logTypes = logTypes
        self.nLogs = {}
        self.Logs  = {}
        self.Info  = {}
        #self.logInfo = ['Process','Filename','Item','Note','Help']

        for t in self.logTypes:
            self.nLogs[t] = 0
            self.Logs[t] = []
           
    #-----------------------------------------------------
    # Adds a standard entry in the Log
    #-----------------------------------------------------
    def add_log(self, logType, logText, logItem, logNote, logHelp):
        lDict = {
            'Process': self.logProcess, 
            'Text': logText, 
            'Item': str(logItem), 
            'Note': logNote, 
            'Help': logHelp,
#            'Time': datetime.now() 
            }
        logType = logType[0].upper()+logType[1:].lower()
        if logType in self.logTypes:
            self.Logs[logType].append(lDict)
            self.nLogs[logType] = self.nLogs[logType] + 1

    #-----------------------------------------------------
    def error(self,logText, logItem, logNote=None, logHelp=None):
        self.add_log('Error',logText,logItem,logNote,logHelp)
    #-----------------------------------------------------
    def warning(self,logText,logItem, logNote=None, logHelp=None):
        self.add_log('Warning',logText,logItem,logNote,logHelp)
    #-----------------------------------------------------
    def info(self,logText,logItem, logNote=None, logHelp=None):
        self.add_log('Info',logText,logItem,logNote,logHelp)


    #-----------------------------------------------------
    # Remove duplicate log enties
    #-----------------------------------------------------
    def select_unique(self,logTypes= LOG_TYPES):
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

    #-----------------------------------------------------
    # Reset the log entries to 0
    #-----------------------------------------------------
    def reset(self):
        self.nLogs = {}
        self.Logs  = {}
        self.Info  = {}

        for t in self.logTypes:
            self.nLogs[t] = 0
            self.Logs[t] = []        

    #-----------------------------------------------------
    def __str__(self):
        return(f"{self.logProcess}")

    #-----------------------------------------------------
    def __repr__(self):
        return(f"{self.logProcess} {self.nLogs}")

    #-----------------------------------------------------
    # Show log entries in logger.info
    #-----------------------------------------------------
    def show(self,logTypes=LOG_TYPES):
        for t in logTypes:
            for l in self.Logs[t]:
                logger.info(f"[{t:7s}] {l['Process']} : {l['Text']} {l['Item']} {l['Note']} () {l['Help']} ")

    #-----------------------------------------------------
    # Show log entries in logger.info
    #-----------------------------------------------------
    def get_nlog(self,logTypes=LOG_TYPES):
        nLog = 0
        for t in logTypes:
            nLog += self.nLogs[t]
        return(nLog)

    #-----------------------------------------------------
    def get_aslist(self,logTypes=LOG_TYPES):
    #-----------------------------------------------------
        retLst = []
        for t in logTypes:
            for l in self.Logs[t]:
                retLst.append({'Type': t, } | l)
        return(retLst)

    #-----------------------------------------------------
    def get_asdf(self,logTypes= LOG_TYPES):
    #-----------------------------------------------------
        return(pd.DataFrame(self.get_aslist(logTypes=logTypes)))

    #-----------------------------------------------------
    def get_ashtml(self,logTypes=LOG_TYPES,classes=None,columns=None,index=False):
    #-----------------------------------------------------
    #     df = self.get_asdf(logTypes=logTypes)
    #     html = df.to_html(columns=columns,classes=classes,index=index).replace("\\n","<br>")
    #     return(html)
        log_data=self.get_asdf(logTypes=logTypes)
        if columns:
            try:
                log_data=log_data[columns]
            except Exception as err:
                raise err
        # Convert the DataFrame's rows to a list of tuples
        table_data = [row for row in log_data.itertuples(index=index)]
        # Convert the DataFrame's columns to a list of strings
        table_header = list(log_data.columns)
        table_dict= {
                    'rows': table_data,
                    'columns': table_header
                    }
        return(table_dict)


    #-----------------------------------------------------
    def info(self,logTypes=LOG_TYPES):
    #-----------------------------------------------------
        self.info={}
        for t in logTypes:
            self.info[t]=[]
            for l in self.Logs[t]:
                note=str(l['note']).replace("'", "").replace('"', '')
                print_info=f"{l['Process']}_{note}_{l['Item']}_{l['Help']}"
                self.info[t].append(print_info) 

    
    #-----------------------------------------------------
    def log_to_UI(self,logTypes=LOG_TYPES):
    #-----------------------------------------------------
        info={} #info=[]
        for t in logTypes:
            # print(f"-- {t.upper():8} ({self.nLogs[t]:3}) ------------------------------------------------------")
            info[t]=[]
            for l in self.Logs[t]:
                print(f"{l['Process']}-{l['Note']} ({l['Item']}) {l['Help']} ")
                description=str(l['Note']).replace("'", "").replace('"', '')
                print_info=f"{l['Process']}_{description}_{l['Item']}_{l['Help']}"
                info[t].append(print_info) # info.append(print_info)
       
        return info

    
    #-----------------------------------------------------
    @classmethod
    def from_aslist(cls, logProcess, logTypes, aslist):
    #-----------------------------------------------------
        instance = cls(logProcess=logProcess, logTypes=logTypes)
        for log in aslist:
            logType = log['Type']
            instance.Logs[logType].append({k: v for k, v in log.items() if k != 'Type'})
            instance.nLogs[logType] += 1
        return instance