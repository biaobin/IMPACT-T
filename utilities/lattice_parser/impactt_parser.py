import sys
import re
from collections import defaultdict
from copy import deepcopy
from math import sqrt, pi
import math
import const
from lattice_parser import lattice_parser
from shutil import copyfile

nest_dict = lambda: defaultdict(nest_dict) # to define a['key1']['key2'] = value
                 

# class: impactt_parser
#======================
class impactt_parser(lattice_parser):
    '''
    transform lte.impt => ImpactT.in.
    
    '''
    def __init__(self, fileName, lineName):
        lattice_parser.__init__(self,fileName, lineName)

        # get brief lte.impz lines
        self.lines = self.get_brieflines()
        
        self.control = {}
        self.beam = {}
        self.lattice = nest_dict()
        
        # initiate with default values
        self.__default_control()
        self.__default_beam()
        self.__default_lattice()
        
        # update with read in *.impz file
        self.control = self.update_control()
        self.beam    = self.update_beam()
        self.lattice = self.update_trackline()

        # some parameters of the beam
        self.Scxl = None
        self.gam0 = None
        self.gambet0 = None
        self.bet0 = None
        self.mass = None
        self.update_paras()
        
    def write_impacttin(self):
        '''
        generate ImpactT.in file.
        '''
        control_file = self.impacttin_control()
        lattice_file = self.impacttin_lattice()
        
        with open('ImpactT.in','w') as f:
            f.write(control_file)
            f.write(lattice_file)

        f.close()     

    # default control, beam, lattice dicts
    #==============================================================================
    def __default_control(self):    
        self.control['CORE_NUM_T'] = 1  
        self.control['CORE_NUM_L'] = 1  
        self.control['DT'] = 1e-12  
        self.control['MAX_STEP'] = 1E6  
        self.control['NBUNCH'] = 1 
        self.control['DIM'] = 6 
        self.control['ERROR'] = 0
        self.control['DIAG'] = 1
        self.control['IMAGE_SC'] = 0
        self.control['IMAGE_STOP_POS'] = 0.02
        self.control['MESHX'] = 32
        self.control['MESHY'] = 32
        self.control['MESHZ'] = 32
        self.control['XRAD'] = 0.015
        self.control['YRAD'] = 0.015
        self.control['PERDLEN'] = 10.0
        self.control['RESTART'] = 0
        self.control['NEMISSION'] = -1
        self.control['TEMISSION'] = 0.0
        self.control['KINETIC_ENERGY'] = 0  # kinetic energy W, E=W+E0
        self.control['BKENERGY'] = None     
        self.control['FREQ_RF_SCALE'] = 2.856e9
        self.control['INI_T'] = 0.0

        # turn all para values to str data type
        for key in self.control:
            self.control[key] = str( self.control[key] )

    def __default_beam(self):
        self.beam['MASS'] = const.electron_mass
        self.beam['CHARGE'] = -1.0
        self.beam['DISTRIBUTION_TYPE'] = 2
        self.beam['NP']   = int(1e3)
        self.beam['TOTAL_CHARGE'] = 1.0e-9 #[C]
        
        self.beam['SIGX'] = 0.0     #[m]
        self.beam['SIGY'] = 0.0     #[m]      
        self.beam['SIGZ'] = 0.0     #[m]
        
        self.beam['SIGPX'] = 0.0    #sig_gambetx/gambet0 [rad]
        self.beam['SIGPY'] = 0.0    #sig_gambety/gambet0 [rad]
        self.beam['SIGPZ'] = 0.0    
        
        # for twiss parameters settings
        self.beam['EMIT_X'] = 0.0
        self.beam['EMIT_NX'] = 0.0   
        self.beam['BETA_X'] = 1.0
        self.beam['ALPHA_X'] = 0.0

        self.beam['EMIT_Y'] = 0.0
        self.beam['EMIT_NY'] = 0.0   
        self.beam['BETA_Y'] = 1.0
        self.beam['ALPHA_Y'] = 0.0

        self.beam['EMIT_Z'] = 0.0
        self.beam['EMIT_NZ']= 0.0   
        self.beam['BETA_Z'] = 1.0
        self.beam['ALPHA_Z']= 0.0

        self.beam['DX']= 0.0
        self.beam['DY']= 0.0
        self.beam['DZ']= 0.0

        # turn all para values to str data type
        for key in self.beam:
            self.beam[key] = str(self.beam[key])

    def __default_lattice(self):
        # drift, 0 element
        self.lattice['DRIFT']['ZEDGE'] = 0.0
        self.lattice['DRIFT']['L'] = 0.0

        # quad, 1 element
        self.lattice['QUAD']['ZEDGE'] = 0.0
        self.lattice['QUAD']['L'] = 0.0
        self.lattice['QUAD']['GRAD'] = 0.0 
        self.lattice['QUAD']['FILEID'] = 0.0 
        self.lattice['QUAD']['DX'] = 0.0
        self.lattice['QUAD']['DY'] = 0.0
        self.lattice['QUAD']['ROTATE_X'] = 0.0 #[rad]
        self.lattice['QUAD']['ROTATE_Y'] = 0.0
        self.lattice['QUAD']['ROTATE_Z'] = 0.0
        self.lattice['QUAD']['FREQ'] = 0.0
        self.lattice['QUAD']['PHASE']= 0.0

        # Sol, 3 element
        self.lattice['SOL']['ZEDGE']=0.0
        self.lattice['SOL']['L']=0.0
        self.lattice['SOL']['FILEID']= None

        # SolRF, 105 element
        self.lattice['SOLRF']['ZEDGE'] = 0.0
        self.lattice['SOLRF']['L']     = 0.0
        self.lattice['SOLRF']['EMAX'] = 1.0
        self.lattice['SOLRF']['FREQ']  = 2856e6
        self.lattice['SOLRF']['PHASE'] = 0.0
        self.lattice['SOLRF']['DX'] = 0.0
        self.lattice['SOLRF']['DY'] = 0.0
        self.lattice['SOLRF']['ROTATE_X'] = 0.0 #[rad]
        self.lattice['SOLRF']['ROTATE_Y'] = 0.0
        self.lattice['SOLRF']['ROTATE_Z'] = 0.0
        self.lattice['SOLRF']['SCALEB']   = 0.0
        self.lattice['SOLRF']['FILEID'] = None
        self.lattice['SOLRF']['Z1'] = None
        self.lattice['SOLRF']['Z2'] = None
        self.lattice['SOLRF']['L_FOURIER_EXP'] = None
        
        # TWS
        self.lattice['TWS']['ZEDGE']=0.0
        self.lattice['TWS']['L']=0.0
        self.lattice['TWS']['CAV_MODE']=2*pi/3
        self.lattice['TWS']['LCOUP']=0.052464
        self.lattice['TWS']['LCAV']=0.104926
        self.lattice['TWS']['EMAX']=0.0
        self.lattice['TWS']['FREQ']=2856E6
        self.lattice['TWS']['PHASE']=0.0
        self.lattice['TWS']['FILEID_1']=None
        self.lattice['TWS']['FILEID_2']=None
        self.lattice['TWS']['FILEID_3']=None
        self.lattice['TWS']['FILEID_4']=None
        self.lattice['TWS']['DX'] = 0.0
        self.lattice['TWS']['DY'] = 0.0
        self.lattice['TWS']['ROTATE_X'] = 0.0 #[rad]
        self.lattice['TWS']['ROTATE_Y'] = 0.0
        self.lattice['TWS']['ROTATE_Z'] = 0.0
        self.lattice['TWS']['SCALEB']   = 0.0
        
        # 112, EMfldCy1
        self.lattice['EMFLDCYL']['ZEDGE'] = 0.0
        self.lattice['EMFLDCYL']['L'] = 0.0
        self.lattice['EMFLDCYL']['FREQ'] = 2856e6
        self.lattice['EMFLDCYL']['PHASE'] = 0.0
        self.lattice['EMFLDCYL']['FILEID'] = None

        # -2 element
        self.lattice['WATCH']['ZEDGE'] = 0.0
        self.lattice['WATCH']['FILENAME_ID'] = 80
        self.lattice['WATCH']['SAMPLE_FREQ'] = 1

        #turn all lattice elem values to string data type
        for elem in self.lattice.keys():
            for key in self.lattice[elem].keys():
                self.lattice[elem][key] = str(self.lattice[elem][key])

    # update the default control, beam, lattice with read-in lte.impz
    #==============================================================================    
    def update_control(self):
        '''
        update control:dict with read in para values.
        '''
        control_sec   = self.get_control_section()
        
        control = deepcopy( self.control )
        for key in control_sec.keys():
            if key in self.control.keys():               
                # update with read in values
                control[key] = control_sec[key]
            else:
                print('Unknown control item:',key,'=',control_sec[key])    
                sys.exit()
        return control
               
    def update_beam(self):
        '''
        update beam:dict with read in para values.
        '''
        beam_sec   = self.get_beam_section()
        
        beam = deepcopy( self.beam )
        for key in beam_sec.keys():
            if key in self.beam.keys():               
                # update with read in values
                beam[key] = beam_sec[key]
            else:
                print('Unknow beam item:',key,'=',beam_sec[key])    
                sys.exit()
        return beam

    def update_trackline(self):
        '''
        update the default track_line lattice para values with lte.impz
        '''
        trackline = self.get_lattice_section() # get the tracking line
        
        j = 0
        for elem in trackline:
            # check if the element type is in self.lattice.keys, i.e. whether in
            # dict_keys(['DRIFT', 'QUAD', 'BEND', 'RFCW', 'WATCH'])          
                         
            # update the not-yet-setting lattice element parameters with the default 
            # values.
            if elem['TYPE'] in self.lattice.keys():
                tmp = deepcopy(self.lattice[elem['TYPE']])
                
                #NAME and TYPR for element in lte.impz are not in the self.lattice[elem['TYPE']].keys()
                table =  list(self.lattice[elem['TYPE']].keys())
                table.append('NAME')
                table.append('TYPE')
                for elem_para in elem.keys():
                    if elem_para not in table:
                        print("ERROR: unknown element parameter",elem_para,"for",elem['NAME'],":",elem['TYPE'])
                        print("PROGRAM STOP!")
                        sys.exit()

                    tmp[elem_para] = elem[elem_para] 
                # replace trackline
                trackline[j] = tmp
                j += 1        
            else:
                print("Unknown element type in lattice section:",elem['TYPE'])
                sys.exit()       
        # turn all elem value to string data type    
        return trackline

    def update_paras(self):  
        Scxl = const.c_light/(2*math.pi*float(self.control['FREQ_RF_SCALE'])) 
        gam0 = (float(self.control['KINETIC_ENERGY'])+float(self.beam['MASS']))/float(self.beam['MASS'])
        gambet0 = sqrt(gam0**2-1.0)
        bet0 = gambet0/gam0

        #update self
        self.Scxl = Scxl
        self.gam0 = gam0
        self.gambet0 = gambet0
        self.bet0 = bet0
        self.mass = float(self.beam['MASS'])
     
    #write in ImpactT.in    
    #==================================================================================================== 
    def impacttin_control(self):
        
        '''
        ImpactZ.in's control section
        '''
        # control section
        #----------------
        control_lines =[]
        #control_lines = [' !=================== \n']
        #control_lines.append('! control section \n')
        #control_lines.append('!=================== \n')
        
        # line-1 
        control_lines.append(self.control['CORE_NUM_T'])
        control_lines.append(self.control['CORE_NUM_L'])
        control_lines.append( '\n' )

        # line-2
        maxstep=int(float(self.control['MAX_STEP']))
        control_lines.append(self.control['DT'])
        control_lines.append(str(maxstep))
        control_lines.append(self.control['NBUNCH'])
        control_lines.append( '\n' )
        
        # line-3
        Np = int(float(self.beam['NP']))
        control_lines.append( self.control['DIM'] )
        control_lines.append( str(Np) )
        control_lines.append( '1' )
        control_lines.append( self.control['ERROR'] )
        control_lines.append( self.control['DIAG'] )
        control_lines.append( self.control['IMAGE_SC'] )
        control_lines.append( self.control['IMAGE_STOP_POS'] )
        control_lines.append( '\n' )
        
        # line-4
        control_lines.append( self.control['MESHX'] )
        control_lines.append( self.control['MESHY'] )
        control_lines.append( self.control['MESHZ'] )
        control_lines.append('1')
        control_lines.append(self.control['XRAD'])
        control_lines.append(self.control['YRAD'])
        control_lines.append(self.control['PERDLEN'])
        control_lines.append('\n')
        
        # line-5
        control_lines.append( self.beam['DISTRIBUTION_TYPE'] )
        control_lines.append( self.control['RESTART'] )
        control_lines.append( '0' )
        control_lines.append( self.control['NEMISSION'] )
        control_lines.append( self.control['TEMISSION'] )
        control_lines.append( '\n' )
        
        # line-6 to line-8
        # beam distribution  
        Scxl    = self.Scxl
        gam0    = self.gam0
        gambet0 = self.gambet0
        bet0    = self.bet0
        
        # (sigi,sigj,sigij) is different from what we know about RMS values,
        # only when sigij=0, sigi refers to RMS size. 
        # set sigij by default 0, for non-zero cases, use twiss-para.
        sigxxp=0.0
        sigyyp=0.0
        sigzzp=0.0

        # If twiss para is given, map twiss to sig
        # ----------------------------------------
        # X-PX
        emitx = float(self.beam['EMIT_X'])
        emit_nx = float(self.beam['EMIT_NX']) 
        if emit_nx != 0.0:
            emitx = emit_nx/gambet0
            
        if emitx != 0.0:
            betax = float(self.beam['BETA_X'])
            alphax = float(self.beam['ALPHA_X'])
            
            # attention: this is not RMS "sigx"
            sigx  = sqrt(emitx*betax/(1+alphax**2))
            sigxp = sqrt(emitx/betax)
            sigxxp= alphax/sqrt(1+alphax**2)
            
            # replace the default values
            self.beam['SIGX'] = str(sigx)
            self.beam['SIGPX'] = str(sigxp)
        
        # Y-PY
        emity = float(self.beam['EMIT_Y'])
        emit_ny = float(self.beam['EMIT_NY']) 
        if emit_ny != 0.0:
            emity = emit_ny/gambet0
            
        if emity != 0.0:           
            betay = float(self.beam['BETA_Y'])
            alphay = float(self.beam['ALPHA_Y'])
            
            sigy   = sqrt(emity*betay/(1+alphay**2))
            sigyp  = sqrt(emity/betay)
            sigyyp = alphay/sqrt(1+alphay**2)
             
            self.beam['SIGY']  = str(sigy)
            self.beam['SIGPY'] = str(sigyp)

        # Z-PZ
        emitz = float(self.beam['EMIT_Z'])
        emit_nz = float(self.beam['EMIT_NZ']) 
        if emit_nz != 0.0:
            emitz = emit_nz/gambet0
            
        if emitz != 0.0:           
            betaz = float(self.beam['BETA_Z'])
            alphaz = float(self.beam['ALPHA_Z'])
            
            sigz   = sqrt(emitz*betaz/(1+alphaz**2))
            sigzp  = sqrt(emitz/betaz)
            sigzzp = alphaz/sqrt(1+alphaz**2)
             
            self.beam['SIGZ']  = str(sigz)
            self.beam['SIGPZ'] = str(sigzp)

        # T-code's z=bet0*ct
        sigz=float(self.beam['SIGZ'])
        dz  =float(self.beam['DZ'])
        self.beam['SIGZ'] = str(bet0*sigz) 
        self.beam['DZ'] = str(bet0*dz)

        # change to IMPACT-T coordinates definition
        #sigX    = float(self.beam['SIGX'])/Scxl  
        #sigY    = float(self.beam['SIGY'])/Scxl  
        #sigZ    = float(self.beam['SIGZ'])/Scxl  
        sigPx   = float(self.beam['SIGPX'])*gambet0
        sigPy   = float(self.beam['SIGPY'])*gambet0        
        sigPz   = float(self.beam['SIGPZ'])*gambet0        
        
        control_lines.append( self.beam['SIGX'] )
        control_lines.append( str(sigPx) )
        control_lines.append( str(sigxxp) )
        control_lines.append( '1.0 1.0' )
        control_lines.append(self.beam['DX'])
        control_lines.append( '0.0 \n' )

        control_lines.append( self.beam['SIGY'] )
        control_lines.append( str(sigPy) )
        control_lines.append( str(sigyyp) )
        control_lines.append( '1.0 1.0' )
        control_lines.append(self.beam['DY'])
        control_lines.append( '0.0 \n' )

        zscale =1.0
        if float(self.beam['DISTRIBUTION_TYPE']) >= 100:
            zscale=1e-9
        
        if self.beam['DISTRIBUTION_TYPE']=='16':
            xmu6 = '0.0'  #xmu6 is used for setting beam energy,
                          #avoid xmu6 being added into the read-in dis.
        else:
            xmu6 = str(gambet0)
        
        control_lines.append( self.beam['SIGZ']    )
        control_lines.append( str(sigPz) ) 
        control_lines.append( str(sigzzp)  )  
        control_lines.append( str(zscale)  )  
        control_lines.append( '1.0' )
        control_lines.append(self.beam['DZ'])
        control_lines.append( xmu6 )
        control_lines.append(' \n')

        # line-9
        current = float(self.beam['TOTAL_CHARGE'])*float(self.control['FREQ_RF_SCALE'])
        current = abs(current)

        if self.control['BKENERGY']=='None':
            self.control['BKENERGY']=self.control['KINETIC_ENERGY']

        control_lines.append(str(current))
        control_lines.append( self.control['BKENERGY'])  
        control_lines.append( self.beam['MASS'] )
        control_lines.append( self.beam['CHARGE'] )
        control_lines.append( self.control['FREQ_RF_SCALE'] )
        control_lines.append( self.control['INI_T'] )
        control_lines.append( '\n' )

        control_lines = ' '.join(control_lines)
        return control_lines

    def impacttin_lattice(self):
        '''
        lattice section in ImpactZ.in file.        
        '''
        lte_lines = []
        #lte_lines = [' !=================== \n']
        #lte_lines.append('! lattice lines \n')
        #lte_lines.append('!=================== \n')        
        
        for elem in self.lattice:
            if elem['TYPE'] == 'DRIFT':
                lte_lines.append( elem['L'] )
                lte_lines.append( '1 1 0' )
                lte_lines.append( elem['ZEDGE'] )
                lte_lines.append('1.0 / \n')

            elif elem['TYPE'] == 'QUAD':
                if elem['FILEID']=='1.0':
                    elem['FILEID']=str(float(elem['L']))
                
                lte_lines.append(elem['L'])
                lte_lines.append('10 20 1')
                lte_lines.append(elem['ZEDGE'])
                lte_lines.append(elem['GRAD'])
                lte_lines.append(elem['FILEID'])
                lte_lines.append('1.0')
                lte_lines.append(elem['DX'])
                lte_lines.append(elem['DY'])
                lte_lines.append(elem['ROTATE_X'])
                lte_lines.append(elem['ROTATE_Y'])
                lte_lines.append(elem['ROTATE_Z'])
                lte_lines.append(elem['FREQ'])
                lte_lines.append(elem['PHASE'])
                lte_lines.append('/ \n')

            elif elem['TYPE'] == 'SOL':
                if elem['FILEID']=='None':
                    print("ERROR: please give the SOL field ID.")
                    sys.exit()
                lte_lines.append(elem['L'])
                lte_lines.append('10 20 3')
                lte_lines.append(elem['ZEDGE'])
                lte_lines.append('0')
                lte_lines.append(elem['FILEID'])
                lte_lines.append('1.0')
                lte_lines.append('0 0 0 0 0 / \n')

            elif elem['TYPE'] == 'SOLRF':
                if elem['FILEID']=='None':
                    print("ERROR: please give the SOLRF field ID.")
                    sys.exit()
                self.update_rfdatax(elem) 

                lte_lines.append(elem['L'])
                lte_lines.append('10 20 105')
                lte_lines.append(elem['ZEDGE'])
                lte_lines.append(elem['EMAX'])
                lte_lines.append(elem['FREQ'])
                lte_lines.append(elem['PHASE'])
                lte_lines.append(elem['FILEID'])
                lte_lines.append('1.01')
                lte_lines.append(elem['DX'])
                lte_lines.append(elem['DY'])
                lte_lines.append(elem['ROTATE_X'])
                lte_lines.append(elem['ROTATE_Y'])
                lte_lines.append(elem['ROTATE_Z'])
                lte_lines.append(elem['SCALEB'])
                lte_lines.append('/ \n')

            elif elem['TYPE'] == 'EMFLDCYL':
                if elem['FILEID']=='None':
                    print("ERROR: please give the EMFLDCYL field ID.")
                    sys.exit()
                lte_lines.append(elem['L'])
                lte_lines.append('10 20 112')
                lte_lines.append(elem['ZEDGE'])
                lte_lines.append('1.01')
                lte_lines.append(elem['FREQ'])
                lte_lines.append(elem['PHASE'])
                lte_lines.append(elem['FILEID'])
                lte_lines.append('1.01 0 0 0 0 0 / \n')


            elif elem['TYPE'] == 'TWS':
                if elem['FILEID_1']=='None':
                    print("ERROR: please give the TWS field ID.")
                    sys.exit()

                x=int(float(elem['FILEID_1']))
                if elem['FILEID_2']=='None':
                    elem['FILEID_2'] = str(x+1)
                if elem['FILEID_3']=='None':
                    elem['FILEID_3'] = str(x+2)
                if elem['FILEID_4']=='None':
                    elem['FILEID_4'] = str(x+3)

                amp=float(elem['EMAX'])
                phi=float(elem['PHASE'])
                dphi=float(elem['CAV_MODE'])
                fac=1/math.sin(dphi)
                dphi1=dphi/pi*180-90
                dphi2=90

                amp12=amp*fac
                phi1 =phi+dphi1
                phi2 =phi+dphi2
               
                tmp=nest_dict()
                tmp['FILEID']=elem['FILEID_1']
                tmp['Z1']='-'+elem['LCOUP']
                tmp['Z2']=elem['LCOUP']
                tmp['L_FOURIER_EXP']=str(2*float(elem['LCOUP']))
                self.update_rfdatax(tmp) 

                tmp['FILEID']=elem['FILEID_2']
                tmp['Z1']='0'
                tmp['Z2']=elem['L']
                tmp['L_FOURIER_EXP']=elem['LCAV']
                self.update_rfdatax(tmp) 

                tmp['FILEID']=elem['FILEID_3']
                tmp['Z1']='0'
                tmp['Z2']=elem['L']
                tmp['L_FOURIER_EXP']=elem['LCAV']
                self.update_rfdatax(tmp) 

                tmp['FILEID']=elem['FILEID_4']
                tmp['Z1']='-'+elem['LCOUP']
                tmp['Z2']=elem['LCOUP']
                tmp['L_FOURIER_EXP']=str(2*float(elem['LCOUP']))
                self.update_rfdatax(tmp) 

                #entrance coupler
                lte_lines.append(elem['LCOUP'])
                lte_lines.append('10 20 105')
                lte_lines.append(elem['ZEDGE'])
                lte_lines.append(elem['EMAX'])
                lte_lines.append(elem['FREQ'])
                lte_lines.append(elem['PHASE'])
                lte_lines.append(elem['FILEID_1'])
                lte_lines.append('1.01')
                lte_lines.append(elem['DX'])
                lte_lines.append(elem['DY'])
                lte_lines.append(elem['ROTATE_X'])
                lte_lines.append(elem['ROTATE_Y'])
                lte_lines.append(elem['ROTATE_Z'])
                lte_lines.append(elem['SCALEB'])
                lte_lines.append('/ \n')

                #cav1
                pos = float(elem['ZEDGE']) + float(elem['LCOUP'])
                lte_lines.append(elem['L'])
                lte_lines.append('10 20 105')
                lte_lines.append(str(pos))
                lte_lines.append(str(amp12))
                lte_lines.append(elem['FREQ'])
                lte_lines.append(str(phi1))
                lte_lines.append(elem['FILEID_2'])
                lte_lines.append('1.01')
                lte_lines.append(elem['DX'])
                lte_lines.append(elem['DY'])
                lte_lines.append(elem['ROTATE_X'])
                lte_lines.append(elem['ROTATE_Y'])
                lte_lines.append(elem['ROTATE_Z'])
                lte_lines.append(elem['SCALEB'])
                lte_lines.append('/ \n')

                #cav1
                lte_lines.append(elem['L'])
                lte_lines.append('10 20 105')
                lte_lines.append(str(pos))
                lte_lines.append(str(amp12))
                lte_lines.append(elem['FREQ'])
                lte_lines.append(str(phi2))
                lte_lines.append(elem['FILEID_3'])
                lte_lines.append('1.01')
                lte_lines.append(elem['DX'])
                lte_lines.append(elem['DY'])
                lte_lines.append(elem['ROTATE_X'])
                lte_lines.append(elem['ROTATE_Y'])
                lte_lines.append(elem['ROTATE_Z'])
                lte_lines.append(elem['SCALEB'])
                lte_lines.append('/ \n')

                #entrance coupler
                pos=pos+float(elem['L'])
                lte_lines.append(elem['LCOUP'])
                lte_lines.append('10 20 105')
                lte_lines.append(str(pos))
                lte_lines.append(elem['EMAX'])
                lte_lines.append(elem['FREQ'])
                lte_lines.append(elem['PHASE'])
                lte_lines.append(elem['FILEID_4'])
                lte_lines.append('1.01')
                lte_lines.append(elem['DX'])
                lte_lines.append(elem['DY'])
                lte_lines.append(elem['ROTATE_X'])
                lte_lines.append(elem['ROTATE_Y'])
                lte_lines.append(elem['ROTATE_Z'])
                lte_lines.append(elem['SCALEB'])
                lte_lines.append('/ \n')
             
            
            elif elem['TYPE']=='WATCH':
                lte_lines.append('0')
                lte_lines.append(elem['SAMPLE_FREQ'])
                lte_lines.append(elem['FILENAME_ID'])
                lte_lines.append('-2')
                lte_lines.append(elem['ZEDGE'])
                lte_lines.append('1.0')
                lte_lines.append(elem['ZEDGE'])
                lte_lines.append('/ \n')

            else:
                print("NOT AVAILABLE ELEMENT TYPE:",elem['TYPE'])
                sys.exit()
       
        lte_lines = ' '.join(lte_lines)
        return lte_lines    
   
    # sub-funcs 
    #===============================================================================
   
    def get_control_section(self):
        '''
        get control section.
        '''
        lines = deepcopy(self.lines)
        
        j=0
        for line in lines:
           lines[j]=line.replace(' ','') 
           j=j+1
           
        # get control section
        pattern1 = re.compile(r'^&control$',re.I)  #not case sensitive
        pattern2 = re.compile(r'^&end$',re.I)
        
        j1,j2 = self._get_index(pattern1, pattern2, lines)
        
        control_lines = lines[j1+1:j2]  # ignored 1st and last element, i.e. &control and &end
        
        control = {}
        for line in control_lines:
            tmp = re.split(';|,',line) # , and ; are both supported
        
            # remove white space
            while '' in tmp:
                tmp.remove('')
            
            for j in tmp:
                tmp2 = j.split('=')
                
                name = tmp2[0].upper()
                value = tmp2[1]
                
                # in case math expression, such as:
                # FREQ_RF_SCALE = 2.998e8/2/pi
                try:
                    eval(value)
                except:
                    pass
                else:
                    value = eval(value.lower())
                    value = str(value)  #back to string
                
                control[name] = value                
                
        return control         
        
    def get_beam_section(self):
        '''
        get beam section.
        '''
        #use deepcopy to keep self.lines unchanged, as white space should be kept in
        #rpn expressions for lattice section  
        lines = deepcopy(self.lines)
        j=0
        for line in lines:
            lines[j] = line.replace(' ','')
            j=j+1
            
        # get lattice section
        pattern1 = re.compile(r'^&beam$',re.I)  #not case sensitive
        pattern2 = re.compile(r'^&end$',re.I)
        
        j1,j2 = self._get_index(pattern1, pattern2, lines)
        
        beam_lines = lines[j1+1:j2]  # ignored 1st and last element, i.e. &beam and &end

        beam = {}
        for line in beam_lines:
            tmp = re.split(';|,',line) # , and ; are both supported
        
            # remove white space
            while '' in tmp:
                tmp.remove('')
            
            for j in tmp:
                tmp2 = j.split('=')
                
                name = tmp2[0].upper()
                value = tmp2[1]
                
                # in case math expression, such as:
                #    total_charge=20*5e-3/2.998e8
                try:
                    eval(value)
                except:
                    pass
                else:
                    value = eval(value.lower())
                    value = str(value)  #back to string
                
                beam[name] = value
                                         
        return beam  

    def get_lattice_section(self):
        '''
        get lattice section. 

        '''   
        lines = self.lines
        
        # get lattice section
        pattern1 = re.compile(r'^&lattice$',re.I)  #not case sensitive
        pattern2 = re.compile(r'^&end$',re.I)
        
        j1,j2 = self._get_index(pattern1, pattern2, lines)
        
        lattice = lines[j1+1:j2]  # ignored 1st and last element, i.e. &lattice and &end
        
        # get the tracked line
        trackline = self.get_trackline(lattice)
        return trackline
               
    def _get_index(self, pattern1, pattern2, lines):
        '''
        get start and finished index of (&control,...,&end)    
        '''       
        j1 = int(1e6)
        j2 = int(1e6)
        cnt = 0
        for line in lines:
            if re.match(pattern1, line):  
                j1 = cnt
                break
            cnt += 1
                    
        cnt = j1
        for j in range(j1, len(lines)):
            if re.match(pattern2,lines[j]):
                j2 = cnt  
                break                  
            cnt += 1    
            
        if j1 == int(1e6):
            print(pattern1,"not found. Input file is wrong.")
            sys.exit()
        elif j2==int(1e6):
            print(pattern2,"not found. Input file is wrong.")
            sys.exit()
                      
        return j1,j2
    
    def update_rfdatax(self, elem):
        '''
        update rfdatax with (z1,z2,L_fourier_exp) 
        '''
        f=open('rfdata'+elem['FILEID'],'r')
        lines=f.readlines()
        f.close()
        if elem['Z1'] != 'None':
            lines[1] = elem['Z1'] + '\n'
        if elem['Z2'] != 'None':
            lines[2] = elem['Z2'] + '\n'
        if elem['L_FOURIER_EXP'] != 'None':
            lines[3] = elem['L_FOURIER_EXP'] + '\n'
            
        if elem['Z1'] !='None' or elem['Z2'] !='None' or elem['L_FOURIER_EXP'] !='None':
            #copyfile('rfdata'+elem['FILEID'],'rfdata'+elem['FILEID']+'.bk')
            
            f=open('rfdata'+elem['FILEID'],'w')
            for line in lines:
                f.write(line)
            f.close()

if __name__=='__main__':
        
    # debug
    # ======
    file_name = 'lte.impt'   
    line_name = 'line'
    
    lte = impactt_parser(file_name,line_name)
    lte.write_impacttin()    
    
        
    
