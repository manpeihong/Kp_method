from tkinter import *
from tkinter import filedialog, messagebox
from tkinter.ttk import Progressbar
import csv
import matplotlib

matplotlib.use("TkAgg") # import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg # implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import numpy as np
import threading
import queue
import cross_platform_config
from sys import platform as _platform
# from scipy.optimize import curve_fit

__version__ = '0.01'


class FIT_FTIR:
    def __init__(self, ks, direction, composition, l_para, m_para, n_para, caltype, listbox,
                 progress_var, wn_beingcalculated):
        self.ks = ks
        self.direction = direction
        self.x = composition
        self.L = l_para
        self.M = m_para
        self.N = n_para

        self.caltype = caltype
        self.listbox = listbox
        self.progress_var = progress_var
        self.wn_beingcalculated = wn_beingcalculated

        self.delta = 0.29

        self.energies = 0
        self.HH = []
        self.LH = []
        self.SO = []
        self.a1 = int(self.direction[0])
        self.a2 = int(self.direction[1])
        self.a3 = int(self.direction[2])
        self.a = np.sqrt(self.a1 * self.a1 + self.a2 * self.a2 + self.a3 * self.a3)
        self.t = 0

        if self.caltype == 1:
            self.calculate()
        elif self.caltype == 0:
            pass

    def cal_energy(self, k):
        energies = []
        kx = k * self.a1 / self.a
        ky = k * self.a2 / self.a
        kz = k * self.a3 / self.a

        H11 = self.L * kx * kx + self.M * (ky * ky + kz * kz)
        H12 = self.N * kx * ky
        H13 = self.N * kx * kz
        H21 = self.N * kx * ky
        H22 = self.L * ky * ky + self.M * (kx * kx + kz * kz)
        H23 = self.N * ky * kz
        H31 = self.N * kx * kz
        H32 = self.N * ky * kz
        H33 = self.L * kz * kz + self.M * (kx * kx + ky * ky)

        def cal_total(Ek):
            _H11 = H11 + k * k + Ek
            _H22 = H22 + k * k + Ek
            _H33 = H33 + k * k + Ek

            total = _H11 * _H22 * _H33 + 2 * H12 * H23 * H13 - _H11 * H23 * H23 - _H22 * H13 * H13 - _H33 * H12 * H12 \
                  - self.delta / 3 * (_H11 * _H22 + _H11 * _H33 + _H22 * _H33 - H12 * H12 - H13 * H13 - H23 * H23)

            return abs(total)

        Ek = 0

        while True:
            if cal_total(Ek + 0.0001) >= cal_total(Ek) and cal_total(Ek - 0.0001) >= cal_total(Ek):
                energies.append(Ek)
                break
            else:
                Ek += 0.0001

        Ek += 0.001

        while True:
            if cal_total(Ek + 0.0001) >= cal_total(Ek) and cal_total(Ek - 0.0001) >= cal_total(Ek):
                energies.append(Ek)
                break
            else:
                Ek += 0.0001

        Ek += 0.2

        while True:
            if cal_total(Ek + 0.0001) >= cal_total(Ek) and cal_total(Ek - 0.0001) >= cal_total(Ek):
                energies.append(Ek)
                break
            else:
                Ek += 0.0001

        return energies

    def calculate(self):
        self.HH = []
        self.LH = []
        self.SO = []

        for i in range(0, len(self.ks)):
            self.progress_var.set(i / len(self.ks) * 100)
            self.wn_beingcalculated.set(self.ks[i])
            k = np.sqrt(self.ks[i]) * 3.7
            result = self.cal_energy(k)
            self.HH.append(-result[0])
            self.LH.append(-result[1])
            self.SO.append(-result[2])

        self.energies = [self.HH, self.LH, self.SO]
        return self.energies

    def returnenergies(self):
        return self.energies

    def addlog(self, string):
        self.listbox.insert(END, string)
        self.listbox.yview(END)


class ThreadedTask1(threading.Thread):
    def __init__(self, queue_1, ks, direction, composition, l_para, m_para, n_para, caltype, listbox,
                 progress_var, wn_beingcalculated):
        threading.Thread.__init__(self)
        self.queue = queue_1
        self.ks = ks
        self.direction = direction
        self.x = composition
        self.L = l_para
        self.M = m_para
        self.N = n_para

        self.caltype = caltype
        self.listbox = listbox
        self.progress_var = progress_var
        self.wn_beingcalculated = wn_beingcalculated

    def run(self):
        fitobject = FIT_FTIR(self.ks, self.direction, self.x, self.L, self.M, self.N, self.caltype, self.listbox,
                             self.progress_var, self.wn_beingcalculated)
        self.queue.put(fitobject.calculate())


class kp_method_GUI(Frame):
    def __init__(self, root, masterroot, listbox, statusbar, status1, status2):
        super().__init__(root, width=cross_platform_config.config.FRAME_WIDTH, bg='#2b2b2b')
        self.root = root
        self.masterroot = masterroot
        self.listbox = listbox
        self.statusbar = statusbar
        self.status1 = status1
        self.status2 = status2
        self.filename = ''
        self.numberofdata = 0
        self.numberofdata2 = 0
        self.colororders = ['blue', 'green', 'cyan', 'magenta', 'yellow', 'black']
        self.colororders2 = ['red', 'green', 'cyan', 'magenta', 'yellow', 'black', 'red', 'green', 'cyan', 'magenta',
                             'yellow', 'black']
        self.wavevector = 0
        self.bandgap = 0
        self.HH1LH1 = 0
        self.LH1SO = 0
        self.Temp = 300
        self.fittype = IntVar()
        self.lowercut = -0.0025
        self.highercut = 0.0025
        self.lowercuty = -1.4
        self.highercuty = 1.4
        self.progress_var = DoubleVar()
        self.text = ''
        self.text2 = ''
        self.wn_beingcalculated = DoubleVar()

        self.calline1 = None
        self.calline2 = None
        self.calline3 = None
        self.energies = []
        self.ks = None
        self.ks_negative = None
        self.k2 = []
        self.energies_0 = []
        self.energies_1 = []
        self.energies_2 = []

        self.osdir = ''

        self.frame0 = Frame(self, width=1000, bg='#262626')
        self.frame0.pack(side=TOP, fill=X, expand=True)
        self.frame0.pack_propagate(0)

        buttonopen = Button(self.frame0, text="Open(⌘+O)",
                            command=self.openfromfile, highlightbackground='#262626', width=9)
        buttonopen.grid(row=0, column=0, columnspan=1)
        buttonsave = Button(self.frame0, text="Save (⌘+S)",
                            command=self.save, highlightbackground='#262626', width=9)
        buttonsave.grid(row=0, column=1, columnspan=1)
        buttonclear = Button(self.frame0, text="Clear(⌘+C)",
                             command=self.clearalldata, highlightbackground='#262626', width=9)
        buttonclear.grid(row=0, column=2, columnspan=1)

        def selectGe():
            self.entry_25.delete(0, END)
            self.entry_25.insert(0, '-32.0')
            self.entry_26.delete(0, END)
            self.entry_26.insert(0, '-5.3')
            self.entry_27.delete(0, END)
            self.entry_27.insert(0, '-32.4')

        Getype = Radiobutton(self.frame0, text="Ge", fg="#a9b7c6", bg='#262626', variable=self.fittype, value=1,
                             command=selectGe)
        Getype.grid(row=0, column=3)
        Getype.select()

        def selectSi():
            self.entry_25.delete(0, END)
            self.entry_25.insert(0, '-7.2')
            self.entry_26.delete(0, END)
            self.entry_26.insert(0, '-3.9')
            self.entry_27.delete(0, END)
            self.entry_27.insert(0, '-7.7')

        Sitype = Radiobutton(self.frame0, text="Si", variable=self.fittype, value=2, command=selectSi, fg="#a9b7c6",
                             bg='#262626')
        Sitype.grid(row=0, column=4)

        self.filepath = Label(self.frame0, text="", bg='#262626', fg="#a9b7c6", width=60)
        self.filepath.grid(row=0, column=5, columnspan=1)

        buttoncal = Button(self.frame0, text="Calculate(⌘+A)",
                           command=self.calculate, highlightbackground='#262626', width=12)
        buttoncal.grid(row=0, column=6, columnspan=1)

        if _platform == "win32" or _platform == "win64":
            buttonopen.config(text='Open(Ct+O)')
            buttonsave.config(text='Save(Ct+S)')
            buttonclear.config(text='Clear(Ct+C)')
            buttoncal.config(text='Calculate(Ct+A)')

        self.frame3 = Frame(self, width=300, bg='#2b2b2b')
        self.frame3.pack(side=RIGHT, fill=BOTH, expand=True)
        self.frame3.pack_propagate(0)

        Label(self.frame3, text='Wavevector:',
              fg="#a9b7c6", bg='#2b2b2b', width=13, anchor=E).grid(row=0, column=0, columnspan=1, sticky=E)
        Label(self.frame3, text='Band gap(eV):',
              fg="#a9b7c6", bg='#2b2b2b', width=13, anchor=E).grid(row=1, column=0, columnspan=1, sticky=E)
        Label(self.frame3, text='HH1-LH1(eV):',
              fg="#a9b7c6", bg='#2b2b2b', width=13, anchor=E).grid(row=2, column=0, columnspan=1, sticky=E)
        Label(self.frame3, text='LH1-SO(ev):',
              fg="#a9b7c6", bg='#2b2b2b', width=13, anchor=E).grid(row=3, column=0, columnspan=1, sticky=E)

        self.wavevector1 = Label(self.frame3, text='{}'.format(self.wavevector), fg="#a9b7c6", bg='#2b2b2b', width=8,
                                 anchor=W)
        self.wavevector1.grid(row=0, column=1, columnspan=1, sticky=E)
        self.bandgap1 = Label(self.frame3, text='{}'.format(self.bandgap), fg="#a9b7c6", bg='#2b2b2b', width=8,
                              anchor=W)
        self.bandgap1.grid(row=1, column=1, columnspan=1, sticky=E)
        self.HH1_LH1 = Label(self.frame3, text='{}'.format(self.HH1LH1), fg="#a9b7c6", bg='#2b2b2b', width=8, anchor=W)
        self.HH1_LH1.grid(row=2, column=1, columnspan=1, sticky=E)
        self.LH1_SO = Label(self.frame3, text='{}'.format(self.LH1SO), fg="#a9b7c6", bg='#2b2b2b', width=8,
                            anchor=W)
        self.LH1_SO.grid(row=3, column=1, columnspan=1, sticky=E)

        label_31 = Label(self.frame3, text='k Low Cut:', width=13, anchor=E, fg="#a9b7c6", bg='#2b2b2b')
        label_32 = Label(self.frame3, text='k High Cut:', width=13, anchor=E, fg="#a9b7c6", bg='#2b2b2b')
        label_33 = Label(self.frame3, text='Energy Low Cut:', width=13, anchor=E, fg="#a9b7c6", bg='#2b2b2b')
        label_332 = Label(self.frame3, text='Energy High Cut:', width=13, anchor=E, fg="#a9b7c6", bg='#2b2b2b')
        self.entry_31 = Entry(self.frame3, highlightbackground='#2b2b2b', width=8)
        self.entry_32 = Entry(self.frame3, highlightbackground='#2b2b2b', width=8)
        self.entry_33 = Entry(self.frame3, highlightbackground='#2b2b2b', width=8)
        self.entry_332 = Entry(self.frame3, highlightbackground='#2b2b2b', width=8)
        label_31.grid(row=5, column=0, sticky=E)
        label_32.grid(row=4, column=0, sticky=E)
        label_33.grid(row=7, column=0, sticky=E)
        label_332.grid(row=6, column=0, sticky=E)
        self.entry_31.grid(row=5, column=1)
        self.entry_32.grid(row=4, column=1)
        self.entry_33.grid(row=7, column=1)
        self.entry_332.grid(row=6, column=1)
        self.entry_31.insert(0, self.lowercut)
        self.entry_32.insert(0, self.highercut)
        self.entry_33.insert(0, self.lowercuty)
        self.entry_332.insert(0, self.highercuty)

        def getbutton31():
            self.entry_31.delete(0, END)
            self.entry_31.insert(0, "%.4f" % self.xclick)

        button31 = Button(self.frame3, text="Get",
                          command=getbutton31, highlightbackground='#2b2b2b', anchor=W, width=3)
        button31.grid(row=4, column=2, sticky=W)

        def getbutton32():
            self.entry_32.delete(0, END)
            self.entry_32.insert(0, "%.4f" % self.xclick)

        button32 = Button(self.frame3, text="Get",
                          command=getbutton32, highlightbackground='#2b2b2b', anchor=W, width=3)
        button32.grid(row=5, column=2, sticky=W)

        def getbutton33():
            self.entry_33.delete(0, END)
            self.entry_33.insert(0, "%.4f" % self.yclick)

        button33 = Button(self.frame3, text="Get",
                          command=getbutton33, highlightbackground='#2b2b2b', anchor=W, width=3)
        button33.grid(row=6, column=2, sticky=W)

        def getbutton332():
            self.entry_332.delete(0, END)
            self.entry_332.insert(0, "%.4f" % self.yclick)

        button332 = Button(self.frame3, text="Get",
                           command=getbutton332, highlightbackground='#2b2b2b', anchor=W, width=3)
        button332.grid(row=7, column=2, sticky=W)

        def CUT():
            self.lowercut = float(self.entry_31.get())
            self.highercut = float(self.entry_32.get())
            self.lowercuty = float(self.entry_33.get())
            self.highercuty = float(self.entry_332.get())
            self.FTIRplot.set_xlim([self.lowercut, self.highercut])
            self.FTIRplot.set_ylim([self.lowercuty, self.highercuty])
            self.canvas.show()

        def zoomall():
            self.entry_31.delete(0, END)
            self.entry_32.delete(0, END)
            self.entry_33.delete(0, END)
            self.entry_332.delete(0, END)
            self.entry_31.insert(0, 400)
            self.entry_32.insert(0, 6000)
            self.entry_33.insert(0, 0)
            self.entry_332.insert(0, 70)
            self.wavenumbers_cut = self.wavenumbers
            self.trans_cut = self.transmissions
            CUT()

        button34 = Button(self.frame3, text="Zoom all",
                          command=zoomall, highlightbackground='#2b2b2b', width=8)
        button34.grid(row=8, column=1)

        button35 = Button(self.frame3, text="CUT",
                          command=CUT, highlightbackground='#2b2b2b', anchor=W, width=3)
        button35.grid(row=8, column=2, sticky=W)

        # seperateline = Label(self.frame3, text='-'*44, anchor=E, fg="#a9b7c6", bg='#2b2b2b')
        # seperateline.grid(row=9, column=0, columnspan=3)

        Label(self.frame3, text='Growth Direction:',
              fg="#a9b7c6", bg='#2b2b2b', width=13, anchor=E).grid(row=9, column=0, sticky=E)
        self.entry_23 = Entry(self.frame3, highlightbackground='#2b2b2b', width=8)
        self.entry_23.grid(row=9, column=1)
        self.entry_23.insert(0, '100')

        Label(self.frame3, text='Composition x:',
              fg="#a9b7c6", bg='#2b2b2b', width=13, anchor=E).grid(row=10, column=0, sticky=E)
        self.entry_24 = Entry(self.frame3, highlightbackground='#2b2b2b', width=8)
        self.entry_24.grid(row=10, column=1)
        self.entry_24.insert(0, '0.21')

        def getbutton22():
            self.entry_24.delete(0, END)
            self.entry_24.insert(0, "%.4f" % self.composition)

        button22 = Button(self.frame3, text="Get",
                          command=getbutton22, highlightbackground='#2b2b2b', anchor=W, width=3)
        button22.grid(row=10, column=2, sticky=W)

        Label(self.frame3, text='L:',
              fg="#a9b7c6", bg='#2b2b2b', width=13, anchor=E).grid(row=11, column=0, sticky=E)
        self.entry_25 = Entry(self.frame3, highlightbackground='#2b2b2b', width=8)
        self.entry_25.grid(row=11, column=1)
        self.entry_25.insert(0, '-32.0')
        self.checkcapd = IntVar()
        checkbox7 = Checkbutton(self.frame3, text="", variable=self.checkcapd, bg='#2b2b2b')
        checkbox7.grid(row=11, column=2, sticky=W)
        checkbox7.select()

        Label(self.frame3, text='M:',
              fg="#a9b7c6", bg='#2b2b2b', width=13, anchor=E).grid(row=12, column=0, sticky=E)
        self.entry_26 = Entry(self.frame3, highlightbackground='#2b2b2b', width=8)
        self.entry_26.grid(row=12, column=1)
        self.entry_26.insert(0, '-5.3')
        self.checkmctcapd = IntVar()
        checkbox8 = Checkbutton(self.frame3, text="", variable=self.checkmctcapd, bg='#2b2b2b')
        checkbox8.grid(row=12, column=2, sticky=W)
        checkbox8.select()

        Label(self.frame3, text='N:',
              fg="#a9b7c6", bg='#2b2b2b', width=13, anchor=E).grid(row=13, column=0, sticky=E)
        self.entry_27 = Entry(self.frame3, highlightbackground='#2b2b2b', width=8)
        self.entry_27.grid(row=13, column=1)
        self.entry_27.insert(0, '-32.4')
        self.checkbased = IntVar()
        checkbox9 = Checkbutton(self.frame3, text="", variable=self.checkbased, bg='#2b2b2b')
        checkbox9.grid(row=13, column=2, sticky=W)
        checkbox9.select()

        self.frame1 = Frame(self, width=700, bg='#2b2b2b')
        self.frame1.pack(side=LEFT, fill=BOTH, expand=True)
        # self.frame1.pack_propagate(0)

        self.FTIRfigure = Figure(figsize=(7, 6), dpi=100)
        self.FTIRfigure.subplots_adjust(left=0.05, bottom=0.05, right=0.92, top=0.95)
        self.FTIRplot = self.FTIRfigure.add_subplot(111)

        self.FTIRplot.plot(self.k2, self.energies_0)
        self.FTIRplot.set_xlim([self.lowercut, self.highercut])
        self.FTIRplot.set_ylim([self.lowercuty, self.highercuty])
        self.FTIRplot.spines['left'].set_position('zero')
        # self.FTIRplot.spines['right'].set_color('none')
        self.FTIRplot.spines['bottom'].set_position('zero')
        # self.FTIRplot.spines['top'].set_color('none')
        self.FTIRplot.xaxis.set_ticks(np.arange(self.lowercut, self.highercut, 0.0005))
        self.FTIRplot.xaxis.set_ticks_position('bottom')
        self.FTIRplot.yaxis.set_ticks(np.arange(self.lowercuty, self.highercuty + 0.2, 0.2))
        self.FTIRplot.yaxis.set_ticks_position('left')
        self.FTIRplot.grid(True)

        self.vline = self.FTIRplot.axvline(x=0, visible=True, color='k', linewidth=0.7)
        self.hline = self.FTIRplot.axhline(y=0, visible=True, color='k', linewidth=0.7)
        self.dot = self.FTIRplot.plot(0, 0, marker='o', color='r')

        self.titleplot = self.FTIRplot.twinx()
        self.titleplot.set_xlabel('k2 (atomic units)')
        self.titleplot.set_ylabel('Energy (meV)')
        self.titleplot.yaxis.set_ticks([])

        # a tk.DrawingArea
        self.canvas = FigureCanvasTkAgg(self.FTIRfigure, self.frame1)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.frame1)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        def on_key_event(event):
            # print('you pressed %s' % event.key)
            key_press_handler(event, self.canvas, self.toolbar)

        self.canvas.mpl_connect('key_press_event', on_key_event)
        self.canvas.mpl_connect('button_press_event', self.onpick)

        def openfromfile_event(event):
            self.openfromfile()

        def calculate_event(event):
            self.calculate()

        def clearalldata_event(event):
            self.clearalldata()

        if _platform == "darwin" or _platform == "linux" or _platform == "linux2":
            masterroot.bind('<Command-o>', openfromfile_event)
            masterroot.bind('<Command-a>', calculate_event)
            masterroot.bind('<Command-c>', clearalldata_event)
        elif _platform == "win32" or _platform == "win64":
            masterroot.bind('<Control-o>', openfromfile_event)
            masterroot.bind('<Control-a>', calculate_event)
            masterroot.bind('<Control-c>', clearalldata_event)

        self.pack()

    def openfromfile(self):
        """Open a band structure .csv file."""
        if self.numberofdata >= 6:
            self.addlog('Cannot add more data.')
            return
        self.wavenumbers = []
        self.transmissions = []
        self.absorptions = []
        file = filedialog.askopenfile(mode='r', defaultextension=".csv")
        if file is None:  # asksaveasfile return `None` if dialog closed with "cancel".
            return
        self.filename = file.name

        i = -1
        while self.filename[i] != '/':
            i -= 1
        self.filename = self.filename[i + 1:None]

        if self.filename[-4:None] != ".csv" and self.filename[-4:None] != ".CSV":
            self.addlog('{} format is not supported. Please select a .CSV file to open.'.format(self.filename[-4:None]))
            return

        self.filepath.config(text=self.filename)
        if self.filename[0:3] == 'Abs':
            with open(file.name, 'r') as f:
                reader = csv.reader(f, delimiter=',')
                for row in reader:
                    try:
                        self.wavenumbers.append(float(row[0]))
                        self.absorptions.append(float(row[1]))
                    except ValueError:
                        pass
            file.close()

            self.absorptionplot = self.FTIRplot.twinx()
            self.fitline_absorption = self.absorptionplot.plot(self.wavenumbers, self.absorptions, self.colororders2[self.numberofdata2], label=self.filename)
            self.absorptionplot.set_ylabel('Absorption Coefficient (cm-1)')
            self.absorptionplot.set_xlim([self.lowercut, self.highercut])
            self.absorptionplot.set_ylim([0, 10000])

            legend = self.absorptionplot.legend(loc='upper right', shadow=True)
            frame = legend.get_frame()
            frame.set_facecolor('0.90')

            # Set the fontsize
            for label in legend.get_texts():
                label.set_fontsize('medium')

            for label in legend.get_lines():
                label.set_linewidth(1.5)

            self.addlog('Added data {} ({})'.format(self.filename, self.colororders2[self.numberofdata2]))
            self.numberofdata2 += 1

        else:
            with open(file.name, 'r') as f:
                reader = csv.reader(f, delimiter=',')
                for row in reader:
                    try:
                        self.wavenumbers.append(float(row[0]))
                        self.transmissions.append(float(row[1]))
                    except ValueError:
                        pass
            file.close()

            # self.FTIRplot = self.FTIRfigure.add_subplot(111)
            self.FTIRplot.plot(self.wavenumbers, self.transmissions, self.colororders[self.numberofdata], label=self.filename)
            self.FTIRplot.set_xlim([self.lowercut, self.highercut])
            self.FTIRplot.set_ylim([self.lowercuty, self.highercuty])
            self.FTIRplot.set_xlabel('Wavenumbers (cm-1)')
            self.FTIRplot.set_ylabel('Transmission (%)')
            self.FTIRplot.grid(True)

            legend = self.FTIRplot.legend(loc='upper right', shadow=True)
            frame = legend.get_frame()
            frame.set_facecolor('0.90')

            # Set the fontsize
            for label in legend.get_texts():
                label.set_fontsize('medium')

            for label in legend.get_lines():
                label.set_linewidth(1.5)

            # plt.savefig("test.png", dpi=300)
            # self.FTIRfigure.show()

            self.addlog('Added data {} ({})'.format(self.filename, self.colororders[self.numberofdata]))
            self.numberofdata += 1

        self.canvas.show()

        if len(self.wavenumbers) == 5810:
            self.addlog('Sample is probably characterized at EPIR.')
        elif len(self.wavenumbers) == 1946:
            self.addlog('Sample is probably characterized at UIC.')

    def save(self):
        """Save the calculated band structure to a csv file. """
        saveornot = messagebox.askquestion(" ", "Save the result as a .csv file?", icon='warning')
        if saveornot == 'yes':
            saveascsv = filedialog.asksaveasfilename(defaultextension='.csv')
            if saveascsv is None:

                return
            if saveascsv[-4:None] != ".csv" and saveascsv[-4:None] != ".CSV":
                self.addlog('Only .csv file can be saved.')
                return
            f = open(saveascsv, "w")

            for i in range(0, len(self.k2)):
                f.write("{0:.6e},{1:.6e},{2:.6e},{3:.6e}\n".format(self.k2[i], self.energies_0[i],
                                                                   self.energies_1[i], self.energies_2[i]))

            f.close()

            self.addlog('Saved the file to: {}'.format(saveascsv))
        else:
            return

    def calculate(self):
        """Calculate the valence band structure using Kane model (1956). """
        self.ks = np.arange(0.00001, float(self.entry_32.get()), 0.00001)
        self.ks_negative = -self.ks

        self.addprogressbar()
        self.text2 = self.status2.cget("text")

        self.addlog('*' * 60)
        self.addlog("Calculation in process. Please wait...")

        self.queue = queue.Queue()
        ThreadedTask1(self.queue,  self.ks, str(self.entry_23.get()),
                      float(self.entry_24.get()), float(self.entry_25.get()), float(self.entry_26.get()),
                      float(self.entry_27.get()), 0, self.listbox, self.progress_var, self.wn_beingcalculated).start()
        self.master.after(100, self.process_queue_calculate)

    def process_queue_calculate(self):
        """Multithread for self.calculate()."""
        try:
            self.trackwavenumber()
            result = self.queue.get(0)
            # Show result of the task if needed
            self.energies = result

            for i in range(0, len(self.ks_negative)):
                self.k2.append(self.ks_negative[len(self.ks_negative)-i-1])
                self.energies_0.append(self.energies[0][len(self.ks_negative)-i-1])
                self.energies_1.append(self.energies[1][len(self.ks_negative)-i-1])
                self.energies_2.append(self.energies[2][len(self.ks_negative)-i-1])

            self.k2 += self.ks.tolist()
            self.energies_0 += self.energies[0]
            self.energies_1 += self.energies[1]
            self.energies_2 += self.energies[2]

            self.calline1 = self.FTIRplot.plot(self.k2, self.energies_0, "-", color=self.colororders[self.numberofdata], label='HH band')
            self.calline2 = self.FTIRplot.plot(self.k2, self.energies_1, "-.", color=self.colororders[self.numberofdata], label='LH band')
            self.calline3 = self.FTIRplot.plot(self.k2, self.energies_2, ":", color=self.colororders[self.numberofdata], label='SO band')

            self.FTIRplot.set_xlim([self.lowercut, self.highercut])
            self.FTIRplot.set_ylim([self.lowercuty, self.highercuty])
            self.FTIRplot.xaxis.set_ticks(np.arange(self.lowercut, self.highercut, 0.0005))
            self.FTIRplot.xaxis.set_ticks_position('bottom')
            self.FTIRplot.yaxis.set_ticks(np.arange(self.lowercuty, self.highercuty + 0.2, 0.2))
            self.FTIRplot.yaxis.set_ticks_position('left')
            self.FTIRplot.grid(True)

            legend = self.FTIRplot.legend(loc='upper right', shadow=True)
            frame = legend.get_frame()
            frame.set_facecolor('0.90')

            # Set the fontsize
            for label in legend.get_texts():
                label.set_fontsize('medium')

            for label in legend.get_lines():
                label.set_linewidth(1.5)

            self.canvas.show()

            self.numberofdata += 1

            self.addlog('Calculation complete!')
            self.removeprogressbar()
            self.removewavenumber()

        except queue.Empty:
            self.after(100, self.process_queue_calculate)

    def onpick(self, event):
        """Show information when mouse click on the plot. """
        self.xclick = event.xdata
        self.yclick = event.ydata
        self.vline.set_xdata(self.xclick)
        self.hline.set_ydata(self.yclick)
        try:
            self.dot.pop(0).remove()
        except IndexError:
            if self.xclick is not None and self.yclick is not None:
                self.dot = self.FTIRplot.plot(self.xclick, self.yclick, marker='x', color='r')
            return
        if self.xclick is not None and self.yclick is not None:
            self.dot = self.FTIRplot.plot(self.xclick, self.yclick, marker='x', color='r')

        self.canvas.draw()

        if self.xclick is not None and self.yclick is not None:
            self.wavevector = self.xclick
            if len(self.energies_0) != 0:
                for i in range(0, len(self.k2)-1):
                    if self.k2[i+1] > self.wavevector >= self.k2[i]:
                        self.bandgap = 0 - self.energies_0[i]
                        self.HH1LH1 = self.energies_0[i] - self.energies_1[i]
                        self.LH1SO = self.energies_1[i] - self.energies_2[i]

        else:
            self.wavevector = 0
            self.bandgap = 0
            self.HH1LH1 = 0
            self.LH1SO = 0

        self.wavevector1.config(text='{0:.4f}'.format(self.wavevector))
        self.bandgap1.config(text='{0:.4f}'.format(self.bandgap))
        self.HH1_LH1.config(text='{0:.4f}'.format(self.HH1LH1))
        self.LH1_SO.config(text='{0:.4f}'.format(self.LH1SO))

    def clearalldata(self):
        """Clear all calculation results."""
        if self.k2 != 0 or self.filepath.cget('text') != '':
            self.FTIRplot.clear()
            self.lowercut = -0.0025
            self.highercut = 0.0025
            self.lowercuty = -1.4
            self.highercuty = 1.4
            self.filepath.config(text='')
            self.filename = ''

            self.calline1 = None
            self.calline2 = None
            self.calline3 = None
            self.energies = []
            self.ks = None
            self.ks_negative = None
            self.k2 = []
            self.energies_0 = []
            self.energies_1 = []
            self.energies_2 = []
            self.FTIRplot.plot(self.k2, self.energies_0)
            self.FTIRplot.set_xlim([self.lowercut, self.highercut])
            self.FTIRplot.set_ylim([self.lowercuty, self.highercuty])
            self.FTIRplot.spines['left'].set_position('zero')
            # self.FTIRplot.spines['right'].set_color('none')
            self.FTIRplot.spines['bottom'].set_position('zero')
            # self.FTIRplot.spines['top'].set_color('none')
            self.FTIRplot.xaxis.set_ticks(np.arange(self.lowercut, self.highercut, 0.0005))
            self.FTIRplot.xaxis.set_ticks_position('bottom')
            self.FTIRplot.yaxis.set_ticks(np.arange(self.lowercuty, self.highercuty + 0.2, 0.2))
            self.FTIRplot.yaxis.set_ticks_position('left')
            self.FTIRplot.grid(True)

            self.vline = self.FTIRplot.axvline(x=0, visible=True, color='k', linewidth=0.7)
            self.hline = self.FTIRplot.axhline(y=0, visible=True, color='k', linewidth=0.7)
            self.dot = self.FTIRplot.plot(0, 0, marker='o', color='r')

            self.titleplot = self.FTIRplot.twinx()
            self.titleplot.set_xlabel('k2 (atomic units)')
            self.titleplot.set_ylabel('Energy (meV)')
            self.titleplot.yaxis.set_ticks([])
            self.canvas.show()

            self.numberofdata = 0
            self.numberofdata2 = 0
            self.addlog('*' * 60)

    def addlog(self, string):
        """Add log to the log frame."""
        self.listbox.insert(END, string)
        self.listbox.yview(END)

    def addprogressbar(self):
        self.text = self.status1.cget("text")
        self.status1.pack_forget()

        self.progressbar = Progressbar(self.statusbar, variable=self.progress_var, maximum=100)
        self.progressbar.pack(side=LEFT, fill=X, expand=1)

    def removeprogressbar(self):
        self.progressbar.pack_forget()

        self.status1 = Label(self.statusbar, text=self.text, fg='#a9b7c6', bg='#2b2b2b', bd=1, relief=RIDGE)
        self.status1.pack(side=LEFT, fill=X, expand=True)
        self.status1.pack_propagate(0)

    def trackwavenumber(self):
        self.status2.config(text='k2 = {:.5f}'.format(self.wn_beingcalculated.get()))

    def removewavenumber(self):
        self.status2.config(text=self.text2)


def main():
    root = Tk()
    w = cross_platform_config.config.FRAME_WIDTH  # width for the Tk root
    h = cross_platform_config.config.FRAME_HEIGHT  # height for the Tk root
    ws = root.winfo_screenwidth()  # width of the screen
    hs = root.winfo_screenheight()  # height of the screen
    x = (ws / 2) - (w / 2)
    y = (hs / 4) - (h / 4)
    root.geometry('%dx%d+%d+%d' % (w, h, x, y))

    root.wm_title("kp method modeling v. {}".format(__version__))
    root.configure(background='#2b2b2b')

    # Status bar #
    statusbar = Frame(root, bg='#2b2b2b', bd=1, relief=RIDGE)

    authorLabel = Label(statusbar, text='© 11,2017 Peihong Man.', fg='#a9b7c6', bg='#2b2b2b', bd=1, relief=RIDGE,
                        padx=4.2, width=21)
    authorLabel.pack(side=LEFT)
    authorLabel.pack_propagate(0)

    status1 = Label(statusbar, text='Welcome to kp method modeling.', fg='#a9b7c6', bg='#2b2b2b', bd=1, relief=RIDGE)
    status1.pack(side=LEFT, fill=X, expand=True)
    status1.pack_propagate(0)
    status2 = Label(statusbar, text='v. {}'.format(__version__), fg='#a9b7c6', bg='#2b2b2b', bd=1, relief=RIDGE,
                    width=21)
    status2.pack(side=RIGHT)

    statusbar.pack(side=BOTTOM, fill=X)

    # Log frame #
    logFrame = Frame(root, height=100, bd=0, highlightthickness=0, bg='white')
    logFrame.pack(side=BOTTOM, fill=X, expand=False)
    logFrame.pack_propagate(0)

    scrollbar = Scrollbar(logFrame, bg='#393c43', highlightbackground='#393c43', troughcolor='#393c43')
    scrollbar.pack(side=RIGHT, fill=Y)

    listbox = Listbox(logFrame, fg='#a9b7c6', bg='#393c43', bd=0, selectbackground='#262626', highlightthickness=0,
                      yscrollcommand=scrollbar.set)
    listbox.pack(side=LEFT, fill=BOTH, expand=True)
    scrollbar.config(command=listbox.yview)

    kp_method_GUI(root, root, listbox, statusbar, status1, status2)

    listbox.delete(0, END)
    listbox.insert(END, '*' * 60)
    listbox.insert(END, 'Welcome to kp method modeling!')
    listbox.insert(END, 'This is the log file.')
    listbox.insert(END, 'Click to copy to the clipboard.')
    listbox.insert(END, '*' * 60)
    listbox.yview(END)

    root.mainloop()


if __name__ == '__main__':
    main()
