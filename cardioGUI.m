function varargout = cardioGUI(varargin)
% CARDIOGUI MATLAB code for cardioGUI.fig
%      CARDIOGUI, by itself, creates a new CARDIOGUI or raises the existing
%      singleton*.
%
%      H = CARDIOGUI returns the handle to a new CARDIOGUI or the handle to
%      the existing singleton*.
%
%      CARDIOGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CARDIOGUI.M with the given input arguments.
%
%      CARDIOGUI('Property','Value',...) creates a new CARDIOGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cardioGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cardioGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cardioGUI

% Last Modified by GUIDE v2.5 02-Sep-2018 19:06:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cardioGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @cardioGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before cardioGUI is made visible.
function cardioGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cardioGUI (see VARARGIN)

%Initialization of all the led signals of the GUI. LED signals explanation
%in line 94.

I = imread('greysquare.jpg');
axes(handles.ECGReadAxes);
imshow(I);
axes(handles.ECGStimulateAxes);
imshow(I);
axes(handles.ECGCreateVectorAxes);
imshow(I);
axes(handles.ECGActivateAxes);
imshow(I);
axes(handles.ECGSendDataAxes);
imshow(I);

set(handles.ECGBitalinoON, 'Visible','off')%Makes invisible the static text that tells if the Bitalino is well connected
set(handles.ECGBitalinoOFF, 'Visible','on')%Makes visible the static text that tells if the Bitalino is not connected

% PREDETERMINED FILTER DESIGN
Apass = 1;           % Passband Ripple (dB)
Astop = 80;          % Stopband Attenuation (dB)    
order = 8;            % order
Fhp= 15/500;
Flp=200/500;
FpassBS= 45/500;
FstopBS= 55/500;
FpassBP=20/500;
FstopBP=150/500;
filtLP = designfilt('lowpassiir', 'FilterOrder', order, ...
                    'PassbandFrequency', Flp, 'StopbandAttenuation', ...
                    Astop, 'PassbandRipple', Apass, 'DesignMethod', ...
                    'ellip');
                
filtHP = designfilt('highpassiir', 'FilterOrder', order, ...
                    'PassbandFrequency', Fhp, 'StopbandAttenuation', ...
                    Astop, 'PassbandRipple', Apass, 'DesignMethod', ...
                    'ellip');
                
filtBS = designfilt('bandstopiir', 'FilterOrder', order, ...
                    'PassbandFrequency1', FpassBS, 'PassbandFrequency2', ...
                    FstopBS, 'PassbandRipple', Apass, ...
                    'StopbandAttenuation', Astop, 'DesignMethod', 'ellip');
                
filtBP = designfilt('bandpassiir', 'FilterOrder', order, ...
                    'PassbandFrequency1', FpassBP, 'PassbandFrequency2', ...
                    FstopBP, 'StopbandAttenuation1', Astop, ...
                    'PassbandRipple', Apass, 'StopbandAttenuation2', ...
                    Astop, 'DesignMethod', 'ellip');                

% Choose default command line output for cardioGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cardioGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cardioGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% LED SIGNALS
%In order to know if the button is pushed and its function called, the
%value of the object is got in a boolean form. In this way, if the Value is
%1, we know will turn the led signal green.

%% ENTERED VALUE CONTROL
%In order to know if the string entered in an edit box is a number, when
%the object is entered, it is turned to double form and if it is not
%possible it will show an error dialog to let the user know that it is
%needed a numeric value.


% --- Executes on selection change in ECGBrand_menu.
function ECGBrand_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ECGBrand_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ECGBrand_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ECGBrand_menu


% --- Executes during object creation, after setting all properties.
function ECGBrand_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ECGBrand_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ECGGetDeviceButton.
function ECGGetDeviceButton_Callback(hObject, eventdata, handles)
% hObject    handle to ECGGetDeviceButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%discover devices
devices = daq.getDevices %Gets the id of the device
%create a session with the discovered device
global s
s = daq.createSession('ni');
assignin('base','s',s)
handles.s = s;
guidata(gcbo,handles);
guidata(hObject,handles);
%Add one analog output
addAnalogOutputChannel(s,'Dev1',0,'Voltage');

% --- Executes on button press in ECGReleaseDeviceButton.
function ECGReleaseDeviceButton_Callback(hObject, eventdata, handles)
% hObject    handle to ECGReleaseDeviceButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s = evalin('base', 's');
release(s)
disp('Device released')


% --- Executes on button press in ECGActivateButton.
function ECGActivateButton_Callback(hObject, eventdata, handles)
% hObject    handle to ECGActivateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ECGActivateButton

%% Led Signals explanation: Line 94

a = get(hObject,'Value');
if a == 1
     I = imread('greensquare.png');
     axes(handles.ECGActivateAxes);
     imshow(I);
else
     J = imread('redsquare.jpg');
     axes(handles.ECGActivateAxes);
     imshow(J);
end

% --- Executes on button press in ECGReadButton.
function ECGReadButton_Callback(hObject, eventdata, handles)
% hObject    handle to ECGReadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ECGReadButton

global filtHP filtLP filtBP filtBS FS  
%The global variables will be used to import data to the bitalino plot of
%the cardioGUI from the filter GUI

% Led Signals explanation: Line 94
a = get(hObject,'Value');

if (a == 1)
    I = imread('greensquare.png');
    axes(handles.ECGReadAxes);
    imshow(I);

    display_size_in_seconds=10; %This variable defines the size in terms of time for the bitalino plot
    sampling_period=1/1000;     %Samples are sent from bitalino at 1000 Samples/s

    b = Bitalino %Creation of the object Bitalino from the file 'Bitalino.m'
    %For Bitalino Revolution, comment the line above and uncomment the line below.
    %b = Bitalino('btspp://201709185182',1000)
    set(handles.ECGBitalinoON, 'Visible','on') %Makes visible the static text that tells if the Bitalino is well connected
    set(handles.ECGBitalinoOFF, 'Visible','off') %Makes invisible the static text that tells if the Bitalino is not connected

    b.startBackground; %Start background acquisition

    axes(handles.ECGBitalinoPlot); %Specification of the axes that have to plot the Bitalino signal in the GUI

    figure_samples=(1:round(display_size_in_seconds/sampling_period)); %Samples per figure
    
    %Initialization of the plot: There are going to be 2 subplots where one
    %is going to correspond to the raw signal and the other to the filtered
    %signal in volts. Both graphics will be plotted respect time.
    
    hAx(1) = subplot(211);
    hLine(1) = line('XData',figure_samples, 'YData',zeros(size(figure_samples)), 'Color','b', 'Parent',hAx(1));
    ylabel('ECG RAW (digital)','fontsize',7);

    hAx(2) = subplot(212);
    hLine(2) = line('XData',figure_samples, 'YData',zeros(size(figure_samples)), 'Color','b', 'Parent',hAx(2));
    ylabel('ECG (mV)','fontsize',7);

    set(hAx, 'Box','on', 'YGrid','on');

    axis(hAx(1),[0,display_size_in_seconds,0,1030]); %The raw signal is in bits so the maximum amplitude that the axes will reach will be 1030 (6 bits of range)
    axis(hAx(2),[0,display_size_in_seconds,-2,2]); %The raw signal will be passed through the transfer fucntion in order to obtain the values in volts (maximum 2/-2)

    input_data=zeros(size(figure_samples)); %Initialization of the variable that will contain the obtained data in a matrix with the size of figure_samples
    buffer_for_input_data=zeros(1,10); %Variable that will contain the temporary storage of the data during the loop in matrix form

    clear ('input_value'); %The variable that collects the data in every loop should be empty in order to rewrite it.

    % OLD BITALINO
    % The following shows which columns correlate to each channel.
    % 1             2   3   4   5   6   7   8   9   10  11
    % packetNumber  I0  I1  I2  I3  A0  A1  A2  A3  A4  A5
    %                               EMG EDA ECG       
    
    % BITALINO REVOLUTION
    % The following shows which columns correlate to each channel.
    % 1             2   3   4   5   6   7   8   9   10  11
    % packetNumber  I0  I1  I2  I3  A0  A1  A2  A3  A4  A5
    %                                   EMG ECG EDA  
   
    Stored_Data_ECG = fopen('Stored_Data_ECG.txt','w'); %Creation of the file that will contain the generated data
    
    while (a == 1) %a is the variable that controls either the toggle button of the led is ON (1 = green) or OFF (0 = red).

        input_value=b.read; %Collection of the values through the read fucntion of 'Bitalino.m'
        %The input value will have matrix form in which each column will
        %correspong to a channel ( Find the correlation of columns with
        %channels in line 231)
        
        while size(input_value)==0
            input_value=b.read; %Keep reading until values are received
        end

        buffer_for_input_data = input_value(:,8)'; %Selection of the channel (ECG -> Column 8) from the input value
         
        input_data=[input_data,buffer_for_input_data]; %New data is updated without deleting the previous one.
        slen = length(buffer_for_input_data);
        input_data=input_data(slen+1:end); %Removes old samples 

        input_data_X=(0:sampling_period:(size(input_data,2)-1)*sampling_period); %Definition of the timeline

        ECG_raw = input_data; 
        
        ECG_volts = (((((ECG_raw./((2.^10)-1))-0.5) .* 3.3)./ 1100) .* 1000); %Application of the transfer function to obtain the value in mV.
        
        
        %Filter application: In every condition, it must be checked if the
        %tick square is marked. A tag is given to the tick square in order 
        %to obtain its value, if it is marked, value = 1, otherwise, 0.

        %If the variable filter is 1, the data follows an extra path to be 
        %filtered. The signal "ECG_volts" is filtered using the function 
        %"filter" where it is specified the filter transfer function 
        %(a=>numerator,b=>denominator) obtained from the filters' GUI
        %and the incoming signal "ECG_volts".
        
        %High pass:
        if (get(handles.ECGEditFilterHPcheck, 'value')== true)
            ECG_volts = filter(filtHP, ECG_volts);
        end
        
        %Low pass:
        if (get(handles.ECGEditFilterLPcheck, 'value')== true)
            ECG_volts = filter(filtLP, ECG_volts);
        end
        
        %Band pass:
        if (get(handles.ECGEditFilterBPcheck, 'value')== true)
            ECG_volts = filter(filtBP, ECG_volts);
        end
        
        %Band stop
        if (get(handles.ECGEditFilterBScheck, 'value')== true)
            ECG_volts = filter(filtBS, ECG_volts);
        end
        
        %Using the initialized plot, the ECG_raw signal is defined to be
        %plotted in the line 1 and the ECG_volts siganal in the line 2.
        
        set(hLine(1), 'YData',ECG_raw);     %draw raw signal
        set(hLine(1), 'XData',input_data_X);  

        set(hLine(2), 'YData',ECG_volts);     %draw filtered signal
        set(hLine(2), 'XData',input_data_X);   

        drawnow % force MATLAB to flush any queued displays

        a = get(hObject,'Value'); %Checking if the value of the toggle button is 
        %still 1, in order to continue or break the loop
        
        %Data Saving: Using the fprintf function, the data is stored in a
        %txt file in two different columns depending on if the signal is
        %raw or filtered.
        fprintf(Stored_Data_ECG,'%6s %12s\n','x',ECG_volts);
        fprintf(Stored_Data_ECG,'%6.2f %12.8f\n',ECG_raw);

    end
    
    %In this part the button red is pushed again so its value is 0.
    
    fclose(Stored_Data_ECG);
        
    msgbox('Data saved succesfully');
    
    b.stopBackground; %Stop data acquisition function from 'Bitalino.m'.
    delete(b) %Deletion of the bitalino object.
end 

set(handles.ECGBitalinoON, 'Visible','off')%Makes invisible the static text that tells if the Bitalino is well connected
set(handles.ECGBitalinoOFF, 'Visible','on')%Makes visible the static text that tells if the Bitalino is not connected

J = imread('redsquare.jpg');
axes(handles.ECGReadAxes);
imshow(J);
    

% --- Executes on button press in ECGCreateVectorButton.
function ECGCreateVectorButton_Callback(hObject, eventdata, handles)
% hObject    handle to ECGCreateVectorButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ECGCreateVectorButton

clc %Clear text in command window
evalin( 'base', 'clear all' )

global Data Fs

a = get(hObject,'Value');
% Led Signals explanation: Line 94

if a == 1 %a is the variable that controls either the toggle button of the led is ON (1 = green) or OFF (0 = red).
    
    %% Data acquisition: 
    % -str2double: The data is obtained in string form and passed to double
    % -get: Obtantion of the value of interest from the GUI, handles. is
    % used to call the variable with its assigned tag.
    
    % Entered value control explanation: Line 99
    
    tpoweron = (str2double(get(handles.ECGEditTime,'String')))*10^-3
    if isnan (tpoweron)
        errordlg('You must enter numeric values','Invalid Input','modal')
    end
    
    Pon_ampPercent = str2double(get(handles.ECGEditAmplitudePower,'String'))
    if isnan (Pon_ampPercent)
        errordlg('You must enter numeric values','Invalid Input','modal')
    end
    
    Data_ampPercent = str2double(get(handles.ECGEditAmplitudeData,'String'))
    if isnan (Data_ampPercent)
        errordlg('You must enter numeric values','Invalid Input','modal')
    end
    
    Address = str2double(get(handles.ECGEditIDNum,'String'))
    if isnan (Address)
        errordlg('You must enter numeric values','Invalid Input','modal')
    end
    
    Stim_freq = str2double(get(handles.ECGEditFrequency,'String'))
    if isnan (Stim_freq)
        errordlg('You must enter numeric values','Invalid Input','modal')
    end
    
    Stim_dur = (str2double(get(handles.ECGEditPulseWidth,'String')))*10^-6;
    if isnan (Stim_dur)
        errordlg('You must enter numeric values','Invalid Input','modal')
    end
    
    Stim_b = str2double(get(handles.ECGEditNumBursts,'String'))
    if isnan (Stim_b)
        errordlg('You must enter numeric values','Invalid Input','modal')
    end
    
    %% CODE
    
    On_amp = 3.3; %1 in Tektronix
    Off_amp = 0; %-1 in tektronix
    timeperbit0 = 75e-6; %25e-6;
    timeperbit1 = 75e-6; %25e-6;
    flagtime = 200e-6;
      
    % Signal properties
    trama1=12; %Signal made up of 12 bits
    trama0=3; % Number of synchronization bits
    Fs=100e3; 
    assignin('base','Fs',Fs)
    
    handles.Fs = Fs;
    guidata(gcbo,handles);
    
    Pon_amp = (Pon_ampPercent/100)*On_amp;
    Data_amp = (Data_ampPercent/100)*On_amp;
    
    % Synchronizing heading
    for i=1:trama0 %Heading made up of '1' (trama0 times)
        Trama(i,1)=1;
    end
    
    % Information
    %Address 
    Address=Address-1;
    inputData=dec2bin(Address);
    nBits = length(inputData); %lenght of input string (number of bits)
    extras=8-nBits;
    j=1; 
    if extras ~=0
        for i=1:extras
            encodedData(i) = 0;
        end
        for i=extras+1:8
            if inputData(j) == '1'
                encodedData(i) = 1; 
            else
                encodedData(i) = 0;
            end
            j=j+1;
        end
    else
        for i=1:8
            if inputData(i) == '1'
                encodedData(i) = 1; 
            else
                encodedData(i) = 0;
            end
        end
    end  
    Address=encodedData;        
    
    A=[Address];

    % Parity bit
    slength=length(A);
    cont=0;
    for i=1:slength
        temp=A(i);
        if temp ==1
            cont=cont+1; %Counts number of '1's
        end
    end

    if cont==0
        temp=0;
    else
        temp=mod(cont,2);
    end
    paritybit=temp;
    A=[A,paritybit];

    Trama=[Trama;A'];

    % Manchester
    slength=length(Trama);
    j=1;

    % Rest of the heading
    for i=1:slength
        if Trama(i,1)==1 %RISING EDGE
            TramaMan(j,1)=0;%TramaMan(j+newindex,1)=0;
            TramaMan(j+1,1)=1;%TramaMan(j+1+newindex,1)=1;
        else  %FALLING EDGE
            TramaMan(j,1)=1;%TramaMan(j+newindex,1)=1;
            TramaMan(j+1,1)=0;%TramaMan(j+1+newindex,1)=0;
        end
        j=j+2;
    end
    
    %assignin('base','TramaMan',TramaMan) %Shows TramaMan in worskpace -for debug

    % Vector creation
    Data = single(1); 
    poweronsamples=Fs*tpoweron;
    for i=1:poweronsamples
        Data(i,1)=Pon_amp;
    end
 
   
    nBits=length(TramaMan);

    sampbit1=floor(Fs*timeperbit1);%Samples per bit
    sampbit0=floor(Fs*timeperbit0);%Samples per bit

    index=poweronsamples+1;
    index=int16(index);
    for i= 1:nBits
        if TramaMan(i,1) ==1
            for j=index:index+sampbit1
                Data(j,1)= Data_amp;
            end
             index=index+sampbit1+1;
        else
            for j=index:index+sampbit0
                Data(j,1)= Off_amp;
            end 
             index=index+sampbit0+1;
        end
    end

    slen=length(Data);
    man=Data;

    %Waits for the uC to read, decode and decide
    % Past: flagtime=800e-6;
    Flagsamps=flagtime*Fs;
    %Flag 2
    for i=slen+1:Flagsamps+slen
        Data(i)=Off_amp;
    end
    slen=length(Data);
    

    time_Data = (single(linspace(0,length(Data)/Fs,length(Data))))';
    axes(handles.ECGCreateVectorPlot); %Specification of the axes in which the plot is desired.
    plot(time_Data,Data);
    assignin('base','Data',Data); %Shows "Data" in the workspace
       
    
     % Stimulation vector
     Stim_period = 1/Stim_freq; %seconds
     Stim_time = Stim_b * Stim_period; %seconds
     Stim_samples = Stim_time * Fs;
     Stim_vector= single(1);
     for i = 1:floor(Stim_dur*Fs)
         Stim_vector(i,1) = On_amp;
     end
     for j=i:floor((Stim_period-Stim_dur)*Fs)
         Stim_vector(j,1) = Off_amp;
     end

     Stim_total = repmat(Stim_vector, Stim_b, 1); %Creates a vector with all the bursts for stimulation
     %assignin('base','Stim_total',Stim_total); %Shows "Stim_total" in the workspace
     DataStim = [Data;Stim_total];
     assignin('base','DataStim',DataStim); %Shows "DataStim" in the workspace
     
    time_DataStim = (single(linspace(0,length(DataStim)/Fs,length(DataStim))))';
    %assignin('base','time_DataStim',time_DataStim); %Shows "time_DataStim" in the workspace
    axes(handles.ECGDownlinkPlot);
    plot(time_DataStim,DataStim);
    
    %Verification
    I = imread('greensquare.png');
    axes(handles.ECGCreateVectorAxes);
    imshow(I);
    
    
else
     J = imread('redsquare.jpg');
     axes(handles.ECGCreateVectorAxes);
     imshow(J);
end


% --- Executes on button press in ECGStimulateButton.
function ECGStimulateButton_Callback(hObject, eventdata, handles)
% hObject    handle to ECGStimulateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = get(hObject,'Value');
%Led signals explanation: Line 94

if a == 1  %a is the variable that controls either the toggle button of the led is ON (1 = green) or OFF (0 = red).
    
    Fs = evalin('base', 'Fs')
    s = evalin('base', 's');
    s.Rate = Fs; %USB6216 supports Rates from 0.1 to 250000.0 scans/sec
    DataStim = evalin('base', 'DataStim');
    queueOutputData(s,DataStim);
    s.startForeground; % s.startBackground();  % s.stop();

     
     I = imread('greensquare.png');
     axes(handles.ECGStimulateAxes);
     imshow(I);

else
     J = imread('redsquare.jpg');
     axes(handles.ECGStimulateAxes);
     imshow(J);
end

% --- Executes on button press in ECGSendDataButton.
function ECGSendDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to ECGSendDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ECGSendDataButton
% Fs = 100000
% 
% s.Rate = Fs; %-10000; %USB6216 supports Rates from 0.1 to 250000.0 scans/sec
% queueOutputData(s,DataStim);
% s.startForeground; % s.startBackground();  % s.stop();

 I = imread('ticksquare.png');
 axes(handles.ECGSendDataAxes);
 imshow(I);


%% --- Frequency Text Box
function ECGEditFrequency_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditFrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ECGEditFrequency as text
%        str2double(get(hObject,'String')) returns contents of ECGEditFrequency as a double

% --- Executes during object creation, after setting all properties.
function ECGEditFrequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ECGEditFrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- PulseWidth Text Box
function ECGEditPulseWidth_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditPulseWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ECGEditPulseWidth as text
%        str2double(get(hObject,'String')) returns contents of ECGEditPulseWidth as a double

% --- Executes during object creation, after setting all properties.
function ECGEditPulseWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ECGEditPulseWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Number of Bursts Text Box
function ECGEditNumBursts_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditNumBursts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ECGEditNumBursts as text
%        str2double(get(hObject,'String')) returns contents of ECGEditNumBursts as a double

% --- Executes during object creation, after setting all properties.
function ECGEditNumBursts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ECGEditNumBursts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- IDNumber Text Box
function ECGEditIDNum_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditIDNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ECGEditIDNum as text
%        str2double(get(hObject,'String')) returns contents of ECGEditIDNum as a double

% --- Executes during object creation, after setting all properties.
function ECGEditIDNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ECGEditIDNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Time Power On Text Box
function ECGEditTime_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ECGEditTime as text
%        str2double(get(hObject,'String')) returns contents of ECGEditTime as a double

% --- Executes during object creation, after setting all properties.
function ECGEditTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ECGEditTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Amplitude Power On Text Box
function ECGEditAmplitudePower_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditAmplitudePower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ECGEditAmplitudePower as text
%        str2double(get(hObject,'String')) returns contents of ECGEditAmplitudePower as a double

% --- Executes during object creation, after setting all properties.
function ECGEditAmplitudePower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ECGEditAmplitudePower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Amplitude Data Text Box
function ECGEditAmplitudeData_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditAmplitudeData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ECGEditAmplitudeData as text
%        str2double(get(hObject,'String')) returns contents of ECGEditAmplitudeData as a double

% --- Executes during object creation, after setting all properties.

function ECGEditAmplitudeData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ECGEditAmplitudeData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Time to Read Text Box
function ECGEditTimeToRead_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditTimeToRead (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ECGEditTimeToRead as text
%        str2double(get(hObject,'String')) returns contents of ECGEditTimeToRead as a double

%% Data acquisition: 
% -str2double: The data is obtained in string form and passed to double
% -get: Obtantion of the value of interest from the GUI, handles. is
% used to call the variable with its assigned tag.

% Entered value control explanation: Line 99

TimeToRead = str2double(get(handles.ECGEditTimeToRead,'String'))
if isnan (TimeToRead)
    errordlg('You must enter numeric values','Invalid Input','modal')
end

% --- Executes during object creation, after setting all properties.
function ECGEditTimeToRead_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ECGEditTimeToRead (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Amplitude in Volts Text Box
function ECGEditAmplitudeVolts_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditAmplitudeVolts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ECGEditAmplitudeVolts as text
%        str2double(get(hObject,'String')) returns contents of ECGEditAmplitudeVolts as a double

%% Data acquisition: 
% -str2double: The data is obtained in string form and passed to double
% -get: Obtantion of the value of interest from the GUI, handles. is
% used to call the variable with its assigned tag.

% Entered value control explanation: Line 99

AVolts = str2double(get(handles.ECGEditAmplitudeVolts,'String'))
if isnan (AVolts)
    errordlg('You must enter numeric values','Invalid Input','modal')
end

% --- Executes during object creation, after setting all properties.
function ECGEditAmplitudeVolts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ECGEditAmplitudeVolts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ECGEditFilterHPbutton.
function ECGEditFilterHPbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditFilterHPbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EDIT_HP

% --- Executes on button press in ECGEditFilterHPcheck.
function ECGEditFilterHPcheck_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditFilterHPcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ECGEditFilterHPcheck


% --- Executes on button press in ECGEditFilterLPbutton.
function ECGEditFilterLPbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditFilterLPbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EDIT_LP

% --- Executes on button press in ECGEditFilterLPcheck.
function ECGEditFilterLPcheck_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditFilterLPcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ECGEditFilterLPcheck


% --- Executes on button press in ECGEditFilterBPbutton.
function ECGEditFilterBPbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditFilterBPbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EDIT_BP

% --- Executes on button press in ECGEditFilterBPcheck.
function ECGEditFilterBPcheck_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditFilterBPcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ECGEditFilterBPcheck


% --- Executes on button press in ECGEditFilterBSbutton.
function ECGEditFilterBSbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditFilterBSbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EDIT_BS

% --- Executes on button press in ECGEditFilterBScheck.
function ECGEditFilterBScheck_Callback(hObject, eventdata, handles)
% hObject    handle to ECGEditFilterBScheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ECGEditFilterBScheck


% --- Executes during object creation, after setting all properties.
function ECGCreateVectorPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ECGCreateVectorPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ECGCreateVectorPlot
