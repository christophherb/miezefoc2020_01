import re
import numpy as np
import itertools
import matplotlib.pyplot as plt
from persistenthomology import find_peaks
from lmfit import Model
class Datafile:
    def __init__(self,datapath):
        '''
        initialized the dataset needs only the corresponding path
        
        Parameters: 
        datapath
        '''
        self.datapath = datapath #the datapath
        self.data_dict = {} #all data is stored here
        self.parameters = [] #the header of the columns
        self.scan = '' # thescan command
        self.palfile = [] #the lines of the pal file
        self.channels = []
        self.read = False
            
            
    def __add__(self,seconddata):
        self.read_all()
        seconddata.read_all()
        try:
            allpoints = sorted(self.allsteps.tolist()+seconddata.allsteps.tolist())
            
            #print(len(self.allsteps),len(seconddata.allsteps))
            allpoints = np.array([np.array(allpoints[i]) for i in range(len(allpoints)) if i == 0 or allpoints[i] != allpoints[i-1]])

            #print(allpoints)
            total_len = len(allpoints)
            totalchannels = sorted(list(set(self.channels + seconddata.channels)),key = lambda x: (x[-1]))[::-1]
            totalchannels = sorted(totalchannels,key = lambda x: x[0])
            
            
            allreshaped = {key:np.zeros((total_len,len(totalchannels)))/0 for (key,val) in self.reshaped_data.items()}#newresahped dict
            alloriginal = {key:np.zeros((total_len*len(totalchannels))) for (key,val) in self.data_dict.items()}      #newfaltdict
            
            
            
            for fillindex,datapoint in enumerate(allpoints):
                if datapoint.tolist() in self.allsteps.tolist():

                    #print('oof',self.allsteps,datapoint)
                    #check if the point is in data set 1
                    dataindex = [i for i,j in enumerate(self.allsteps) if np.array_equal(j,datapoint)][0]
                    
                    
                    
                    
                    
                    for channelindex,channel in enumerate(totalchannels):
                        if channel in self.channels:
                            #SELF.POL HERE
                            for key in allreshaped:
                                if (key == 'CNTS' or key == 'M1' or key == 'M2') and not np.isnan(allreshaped[key][fillindex,channelindex]):
                                    alloriginal[key][fillindex*len(totalchannels)+channelindex] +=  self.reshaped_data[key][dataindex,self.pol2pal(channel)]
                                    allreshaped[key][fillindex,channelindex] += self.reshaped_data[key][dataindex,self.pol2pal(channel)]
                                else:
                                    allreshaped[key][fillindex,channelindex] = self.reshaped_data[key][dataindex,self.pol2pal(channel)]
                                    alloriginal[key][fillindex*len(totalchannels)+channelindex] =  self.reshaped_data[key][dataindex,self.pol2pal(channel)]
                if datapoint.tolist() in seconddata.allsteps.tolist():
                    dataindex = [i for i,j in enumerate(seconddata.allsteps) if np.array_equal(j,datapoint)][0]
                    for channelindex,channel in enumerate(totalchannels):
                        if channel in seconddata.channels:
                            for key in allreshaped:
                                if (key == 'CNTS' or key == 'M1' or key == 'M2') and not np.isnan(allreshaped[key][fillindex,channelindex]):
                                    allreshaped[key][fillindex,channelindex] += seconddata.reshaped_data[key][dataindex][seconddata.pol2pal(channel)]
                                    alloriginal[key][fillindex*len(totalchannels)+channelindex] +=  seconddata.reshaped_data[key][dataindex,seconddata.pol2pal(channel)]
                                else:
                                    allreshaped[key][fillindex,channelindex] = seconddata.reshaped_data[key][dataindex][seconddata.pol2pal(channel)]
                                    alloriginal[key][fillindex*len(totalchannels)+channelindex] =  seconddata.reshaped_data[key][dataindex,seconddata.pol2pal(channel)]

            self.allsteps = allpoints
            self.channels = totalchannels           
            self.reshaped_data = allreshaped
            self.data_dict = alloriginal
            return self
        except AttributeError:
            print('initdata missing')
            return seconddata
    def __radd__(self,other):
        self.read_all()
        return self
############################################### READ EVERYTHING ##########################################################################
    def read_all(self):
        if not self.read:
            self.read_data()
            self.read_scan()
            self.read_pols()
            self.read = True
############################################### READ THE DATA INTO A DICT WHICH IS LABLLED ACCORDINGLY Ä#####################################
    def read_data(self):
        data_lines = 0 # where numbers are
        meta_dict = {} 
        data_dict = {} 
        scan_params = [] #scan parameters i.e. CNTS , QH/K/L
        units = False
        with open(self.datapath) as f:#the file is opened
            content = f.readlines()
            for line_in,line in enumerate(content):#content is read line per line
                if units == True:
                    scan_params = re.findall('(?<=\s)\w+(?=\s)',line)#all column headers are saved in scan_params
                    metadata = False
                    units = False
                #print(re.search('DATA',line))
                if re.search('DATA_*',line):#once the DATA file is found in the ILL data structure
                    data_lines = line_in#the line index is saved to indidcate that number are coming from now
                    units = True # units is set to True since the next line tells us what the columns mean
                if re.search('COMND:\s*',line):#once the DATA file is found in the ILL data structure
                    self.scan = re.findall('(?<=COMND:\s).*',line)[0]#the line index is saved to indidcate that number are coming from now
                    units = True 
                if re.search('POLAN:\s*',line):
                    self.palfile += [re.findall('(?<=POLAN:).*',line)[0] ]
                if 'PARAM:' in line or 'VARIA:' in line or 'ZEROS:' in line:
                    listpars = line.split(':')[-1].split(',')
                    for i in listpars:
                        param,value =i.split('=')[0].strip(),float(i.split('=')[-1])
                        meta_dict[param] = value
        data = np.genfromtxt(self.datapath, skip_header=data_lines+2)#data is read in a matrix, skip header is told where to start with data
        
        #handling the odd case of only one point
        if (len(data.shape)) == 1:
            data = data.reshape((1,len(data)))
        #populate data dict for easier post processing
        for i,j in enumerate(scan_params):
            data_dict[j] = data[:,i]
        if 'PAL' in data_dict:
            numberchannels = int(np.max(data_dict['PAL']))
            numberpoints = int(np.max(data_dict['PNT']))
            #print(numberchannels)
            reshaped_data = {key:np.zeros((numberpoints,numberchannels))/0 for (key,val) in data_dict.items()}
            for index in range(len(data_dict['PAL'])):
                
                i = int(data_dict['PNT'][index])-1
                j = int(data_dict['PAL'][index])-1
                #print(i,j)
                for key in data_dict:
                    reshaped_data[key][i,j] = data_dict[key][index]
            self.reshaped_data = reshaped_data
        else:
            numberchannels = 1
            numberpoints = int(np.max(data_dict['PNT']))
            reshaped_data = {key:np.zeros((numberpoints,numberchannels))/0 for (key,val) in data_dict.items()}
            for key in data_dict:
                reshaped_data[key][:,0] = data_dict[key][:]
            self.reshaped_data = reshaped_data
        #params are saved and data dict is handed to object
        self.parameters = scan_params
        self.data_dict = data_dict
        self.read_scan()
        self.meta_dict = meta_dict

    def read_scan(self):
        '''
        reads the scan resulting in the file and returns the suspected x and y axis
        
        Parameters:
        
        none
        '''
        capstring = self.scan.upper()
        if 'BS' in capstring:
            whichscan = 'bs'
        if 'SC' in capstring:
            whichscan = 'sc'
        stringlist = capstring.split()
        xaxis = stringlist[1]
        initvektor,stepvektor  = [],[]
        beforeD = True
        for val in stringlist[2:]:
            if 'D' in val:
                beforeD = False
            if beforeD:
                initvektor += [float(val)]
            else:
                if 'TI' in val:
                    break
                try:
                    stepvektor += [float(val)]
                except ValueError:
                    pass
                    
        initvektor,stepvektor,numbersteps = np.array(initvektor),np.array(stepvektor)[:-1],int(stepvektor[-1])
        #print(initvektor,stepvektor,numbersteps)
        self.initpoint = initvektor
        self.allsteps = np.array([i*stepvektor+initvektor for i in range(numbersteps) ])[:int(max(self.data_dict['PNT']))]
        
        if xaxis == 'QH':
            
            possiblesteps = ['QH','QK','QL','EN']
            realsteps = [possiblesteps[i] for i,val in enumerate(stepvektor) if float(val) != 0]
            if not realsteps: #default value if this does not work ie all are 0 
                realsteps = ['QH']
            self.number_steps = numbersteps
            self.guess_x = realsteps[0]
            self.guess_y = 'CNTS'
            
            
            
        if xaxis == 'A3':
            self.guess_x = 'A3'
            self.guess_y = 'CNTS'
            self.number_steps = numbersteps
        
    def read_pols(self):
        '''
        reads in the pal file to determine the polarization of the pal channels
        sets self.channels to [] if no pal file exists
        '''
        channels = []
        for line in self.palfile:
            if 'HX' in line:
                currents = [abs(float(k)) for k in re.findall('(?<=HX\s).*',line)[0].split()]
                direction = ['x','y','z'][currents.index(max(currents))]
            if 'ON' in line:
                sign = '-'
            if  'OF' in line:
                sign = '+'
            if 'CO' in line:
                channels += [direction+sign]
        self.channels = channels
        
    def pol2pal(self,polarization):
        if not self.channels:
            print('no can do')
        elif polarization not in self.channels:
            print('please use one of x+,x-,y+,y-,z+,z-')
            raise 'StupidError'
        else:
            return self.channels.index(polarization)
    #get counts
    def get_values_pol(self,polarization,values='CNTS'):# = self.guess_y):
        if polarization in self.channels:
            index = self.channels.index(polarization)
            return self.reshaped_data[values][:,index]
        
        
        
    def get_nominal_x(self):
        pass
    
    
    def plotstandard(self,scaling = 1,plotout = True):
        xaxis = self.guess_x
        yaxis = self.guess_y
        polarization = ['x-','y-','z-']
        
        xplot = self.reshaped_data[xaxis][:,self.pol2pal('x-')]
        
        xxbar = self.reshaped_data[yaxis][:,self.pol2pal('x-')]
        yybar = self.reshaped_data[yaxis][:,self.pol2pal('y-')]
        zzbar = self.reshaped_data[yaxis][:,self.pol2pal('z-')]
        
        xxmon = self.reshaped_data['M1'][:,self.pol2pal('x-')]
        yymon = self.reshaped_data['M1'][:,self.pol2pal('y-')]
        zzmon = self.reshaped_data['M1'][:,self.pol2pal('z-')]
        
        normxx,normyy,normzz = xxbar/xxmon,yybar/yymon,zzbar/zzmon
         
        errnxx, errnyy, errnzz = xxbar**0.5/xxmon,yybar**0.5/yymon,zzbar**0.5/zzmon
        
        mperp = 2*normxx -normyy-normzz
        myy = normxx-normyy
        mzz = normxx-normzz
        
        mperperr = (4*errnxx**2 + errnyy**2 + errnzz**2)**0.5
        myyerr = (errnxx**2 + errnyy**2)**0.5
        mzzerr = (errnxx**2 + errnzz**2)**0.5
        
        s = scaling
        if plotout:
            fig,ax = plt.subplots(1,2,figsize = (14,7))



            ax[0].errorbar(xplot,normxx*s,errnxx*s,label = r'x$\bar{x}$',linestyle = '-',markeredgewidth = 0)
            ax[0].errorbar(xplot,normyy*s,errnyy*s,label = r'y$\bar{y}$',linestyle = '-',markeredgewidth = 0)
            ax[0].errorbar(xplot,normzz*s,errnzz*s,label = r'z$\bar{z}$',linestyle = '-',markeredgewidth = 0)


            ax[0].set_xlabel(xaxis)
            ax[0].set_ylabel(yaxis)
            ax[0].legend()
            #ax[0].text(0.1,0.5,'$\Delta$E = {} meV'.format(np.median(dat1.reshaped_data['EN'][0])),transform=ax[0].transAxes)
            ax[0].set_title(self.scan)


            ax[1].errorbar(xplot,mperp*s,mperperr*s,label = r'$M\perp$',linestyle = '-',markeredgewidth = 0)
            ax[1].errorbar(xplot,myy*s,myyerr*s,label = r'$M_{yy}$',linestyle = '-',markeredgewidth = 0)
            ax[1].errorbar(xplot,mzz*s,mzzerr*s,label = r'$M_{zz}$',linestyle = '-',markeredgewidth = 0)


            ax[1].set_xlabel(xaxis)
            ax[1].set_ylabel(yaxis)
            ax[1].legend()
            #ax[0].text(0.1,0.5,'$\Delta$E = {} meV'.format(np.median(dat1.reshaped_data['EN'][0])),transform=ax[0].transAxes)
            ax[1].set_title(self.scan)
        return  xplot,mperp*s,myy*s,mzz*s, mperperr*s,myyerr*s, mzzerr*s 
    
    def fit_gauss(self,xvals,yvals,errs,xlabel='huh',ylabel = 'CNTS',plotout = True, linearoffset = False, 
                  xmin=False,xmax = False, peaks = None,vlines = True,init_fit = False):
        '''
        Fits a set number of gaussians to a curve
        
        Paramerters:
        
        xaxis: Name of the variable to appear on the xaxis
        yaxis: Name of the variable to appear on the yaxis usually counts
        xmin:  Minimum value of x 
        xmax:  Maximum value of x
        peaks: Either list of expected xvalues of peaks or the number of peaks to be found
        vlines: If True, the peakpositions are highlighted in the resulting image
        init_fit: If True, the initial guess is shown, mostly to see how good the initial guess is 
        '''
        def offset(x,c,m = 0):
            '''
            offsetfunction with linear slope returning m*x +c
            
            Parameters:
            
            x: xvalue
            c: constant offset
            m: slope
            '''
            return c + m*x
        def gauss(x,area,center,sigma):
            '''
            normalized gaussian 
            
            Parameters:
            
            x: xvalue
            area: area under the gaussian cruve also height/sigma/(2pi)**0.5
            center: horizontal offset of gaussian
            sigma: well, sigma
            c: unused vertical offset
            
            '''
            return area* np.exp(-(x-center)**2/(2*sigma**2))/(sigma*(2*np.pi)**0.5)
        
        
        x = xvals
        y = yvals
        # x and y value extracted from the data files, necessary to plot the data later
        y_err = errs
        # data to be fitted, will be cropped as necessary
        x_fit = x
        y_fit = y
        y_err_fit = y_err
        #if normalized == True:
        #    y = self.data_dict['NORM']
        #    y_err = self.data_dict['NORMERR']
        #    y_fit = y
        #    y_err_fit = y_err
        # fitting data is reduced to the range from xmin - xmax 
        if np.logical_and(xmin,xmax): 
            x_fit = x[np.where(np.logical_and(x >= xmin,x <= xmax))]
            y_fit = y[np.where(np.logical_and(x >= xmin,x <= xmax))]
            y_err_fit = y_err[np.where(np.logical_and(x >= xmin,x <= xmax))]
            #print(x_fit.shape, y_fit.shape, y_err_fit.shape)
            
        # data to plot a smooth gaussian from xmin to xmax 
        if plotout:
            x_plot = np.linspace(np.min(x_fit),np.max(x_fit),1000)
            
        #finding the x values of the most important peaks if only the number of peaks is given
        if type(peaks) is not list:
            self.peaks,self.peaky = find_peaks(peaks,x_fit,y_fit)
            #print(self.peaks)
        #if a list of xvalues is provided, find the corresponding y values, for peak height determination 
        else:
            self.peaks = peaks
            self.peaky = np.zeros(len(self.peaks))
            
            for i,j in enumerate(self.peaky):
                self.peaky[i] = y_fit[np.where((x_fit- self.peaks[i])**2 == np.min((x_fit- self.peaks[i])**2))][0]
        
        
        if True> 0:#this I dont really know, what is it's purpose, what are its wishes
            paramsi = False
            final = Model(offset,prefix = 'o0_')#final Model, to which all the gaussians are added
            params = final.make_params() 
            params.add('o0_c', value= np.min(y_fit),min = -100)
            
            #if linear offset is enabled, the slope of the offsetfunction is variable, else ti es set to zeor
            if linearoffset:
                params.add('o0_m', value = 0)
            else:
                params.add('o0_m', value = 0, vary = False)
            #print(x_fit,x_fit[9]-x_fit[8])
            for i,peak in enumerate(self.peaks):
                #one gaussian Model is added per peak
                final += Model(gauss,prefix = 'g{}_'.format(i))
                #sigma determination is shit, I will look up a better algorhithm maybe quick integration and hight?
                sigma = abs(x_fit[1]-x_fit[0])
                area = (self.peaky[i]-np.min(y_fit))*(sigma*(2*np.pi)**0.5)#same problem as above
                #print(area,sigma,self.peaks[i])
                #all gaussian params are added as name, value, vary, min, max, expression
                params.add_many(('g{}_area'.format(i), area  ,  True, 0, None,  None),\
                                 ('g{}_center'.format(i),   self.peaks[i],  True,  None, None,  None),\
                                   ('g{}_sigma'.format(i),   sigma,  True,  0.01, 22,  None),\
                               )
            #print(x_fit,y_fit,y_err_fit)
            #Model is fitted to the data with 1/err weights
            res = final.fit(x = x_fit,data = y_fit, params = params, weights = 1/y_err_fit)
            #print('fuck',res.params['g0_center'].stderr,res.params['g1_center'].stderr)
            #res.params.pretty_print()
            #if wanted all fits are printed out
            if plotout:
                
                #plots flat x and y
                self.ax.errorbar(x,y,y_err)
                
                #can plot the initial guess
                if init_fit:
                    self.ax.plot(x_fit, res.init_fit, 'k--')
                
                #stores the plotvals of the individual gaussians in a dict
                plot_result= {}
                #populate the dict and plot the individual gaussians
                for i,_ in enumerate(self.peaks):
                    plot_result[i] = gauss(x_plot, res.params['g{}_area'.format(i)].value, res.params['g{}_center'.format(i)].value, res.params['g{}_sigma'.format(i)].value)
                    self.ax.plot(x_plot,plot_result[i],'b-')
                
                #plot the offset function
                self.ax.plot(x_plot, offset(x_plot,res.params['o0_c'].value,res.params['o0_m']),'g-')
                
                #the complete fit
                self.ax.plot(x_plot,offset(x_plot,res.params['o0_c'].value,res.params['o0_m']) +np.sum([plot_result[i] for i in range(len(self.peaks))],axis = 0),'r--',\
                             label = 'Fit with red. chisqr {:.4}'.format(res.redchi))
                
                #plot vertical lines at the peakpositions if vlines = True
                if vlines:
                    for i,peak in enumerate(self.peaks):
                        center,err,sigma = res.params['g{}_center'.format(i)].value,res.params['g{}_center'.format(i)].stderr,\
                        res.params['g{}_sigma'.format(i)].value 
                        plt.axvline(res.params['g{}_center'.format(i)].value,linestyle = '-',\
                        color = 'black',marker = 'None',linewidth = 1,dashes =[4,4])
                        self.ax.text(res.params['g{}_center'.format(i)].value+0.05*(np.max(x_fit)-np.min(x_fit)), 0.8*np.max(y_fit),\
                        'center ={:.4} \n err = {:.4} \n sigma= {:.4}'.format(center,err,sigma))
                
                #standard stuff, legend, labels
                #print(res.params['g{}_center'.format(i)].stderr)
                self.ax.legend()
                self.ax.set_ylabel(ylabel)
                self.ax.set_xlabel(xlabel)
                plt.tight_layout()
                plt.show()
                #returns the area of the main gausspeak
        return res