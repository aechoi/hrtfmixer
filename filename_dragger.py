import wx

# idk what I'm doing, all I can say is that I'm resting easy
# knowing that one day I will potentially just be redoing this
# whole 'GUI' thing with a single package one day. This is 
# very just cobbled together and is almost certainly bad practice.

# but, i mean, it sorta works for now

app = wx.App(False)


########################################################################
class MyFileDropTarget(wx.FileDropTarget):
    """"""

    #----------------------------------------------------------------------
    def __init__(self, window):
        """Constructor"""
        wx.FileDropTarget.__init__(self)
        self.window = window
        self.fileNames = []

    #----------------------------------------------------------------------
    def OnDropFiles(self, x, y, filenames):
        """
        When files are dropped, write where they were dropped and then
        the file paths themselves
        """
        self.window.SetInsertionPointEnd()
        # self.window.updateText("\n%d file(s) dropped at %d,%d:\n" %
        #                       (len(filenames), x, y))
        for filepath in filenames:
            self.window.updateText('Loaded ' + filepath.split('\\')[-1] + '\n')    
        # print (filenames)
        # return filenames
        self.fileNames = filenames
        frame.Close()
        return True

########################################################################
class DnDPanel(wx.Panel):
    """"""

    #----------------------------------------------------------------------
    def __init__(self, parent):
        """Constructor"""
        wx.Panel.__init__(self, parent=parent)

        self.file_drop_target = MyFileDropTarget(self)
        lbl = wx.StaticText(self, label="Drag a 16-bit .wav file here.\nThe modified file will be saved in the same location.")
        self.fileTextCtrl = wx.TextCtrl(self,
                                        style=wx.TE_MULTILINE|wx.HSCROLL|wx.TE_READONLY)
        self.fileTextCtrl.SetDropTarget(self.file_drop_target)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(lbl, 0, wx.ALL, 5)
        sizer.Add(self.fileTextCtrl, 1, wx.EXPAND|wx.ALL, 5)
        self.SetSizer(sizer)

    #----------------------------------------------------------------------
    def SetInsertionPointEnd(self):
        """
        Put insertion point at end of text control to prevent overwriting
        """
        self.fileTextCtrl.SetInsertionPointEnd()

    #----------------------------------------------------------------------
    def updateText(self, text):
        """
        Write text to the text control
        """
        self.fileTextCtrl.WriteText(text)

########################################################################
class DnDFrame(wx.Frame):
    """"""

    #----------------------------------------------------------------------
    def __init__(self):
        """Constructor"""
        wx.Frame.__init__(self, parent=None, title='Select a File')
        self.panel = DnDPanel(self)
        self.Show()

#########################################################################
frame = DnDFrame()
class dragStarter():
    def __init__(self):
        app.MainLoop()
        self.fileName = frame.panel.file_drop_target.fileNames
#----------------------------------------------------------------------
if __name__.endswith("__main__"):
    # app = wx.App(False)
    # frame = DnDFrame()
    app.MainLoop()
    print(frame.panel.file_drop_target.fileNames)