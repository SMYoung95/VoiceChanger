classdef voiceChangerPlugin < audioPlugin
    
    properties
        % Name of the sliders
        
        % Ring Modulation
        % Frequency of the synthesized wave
        rModFreq = 400;
        % Phase of the synthesized wave
        waveOffset = 0;
        % Wave Type
        waveType = 'sine';
        % Counter to for oscillator to know where to start.
        k = 1;
        % General Volume
        volume = 100;
        Center = 0.05;
        Range = 100;        
        
        % Delay
        % Duration of delay in ms 
        fbDelayNoteDuration = 0;
        % Volume of delay
        WetDryMix = 0;        

        % Filters
        % High and Lowpass types (lowpass(0) , highpass(1))
        filtType = 'Low';
        
        % Pass Frequency of filters (cut-OffFrequency)
        Fc = 3000;
        % Gain of filters
        Q = 6;
        
        finalFreq = 0
        

        
        % Not sure, gotta check up, something to do with the filters
        z = zeros(2);
        
        pitchChange = 0;
       
        
    end
    
    properties (Constant)
        % Initializing the plugin slider, parameters has to match the name
        % of the slider, the relationship between the property and the 
        % position of the GUI control on the plugin's dialog, {'linear'
        % ranges the mapping goes between.
        PluginInterface = audioPluginInterface( ...
            'InputChannels', 2,...
            'OutputChannels', 2,...
            'PluginName', 'Youngs Voice Changer',...
            audioPluginParameter(...
            'rModFreq', 'DisplayName', 'Frequency', 'Mapping', {'lin', 1,1000}), ...
            audioPluginParameter(...
            'waveType', 'DisplayName', 'Wave Type', 'Mapping', {'enum', ...
            'sine', 'sawtooth', 'square'}), ...    
            audioPluginParameter(...
            'volume', 'DisplayName', 'Volume', 'Mapping', {'lin', 0,100}),...
            audioPluginParameter(...
            'fbDelayNoteDuration', 'DisplayName', 'Feedback Delay Echo Duration',...
            'Mapping', {'lin', 0,10000}),...
            audioPluginParameter(...
            'WetDryMix', 'DisplayName', 'WetDryMix',...
            'Mapping', {'lin', 0,10}),...
            audioPluginParameter(...
            'filtType', 'DisplayName', 'Low or HighPass Filter',...
            'Mapping', {'enum', 'Low', 'High'}),...
            audioPluginParameter(...
            'Fc', 'DisplayName', 'Cutoff Freq',...
            'Label', 'Hz', 'Mapping', {'log', 50,5000}),...
            audioPluginParameter(...
            'Q', 'Mapping', {'lin', 0.1,6}),...
            audioPluginParameter(...
            'Range', 'Mapping', {'lin', 0.1,100}),...
            audioPluginParameter(...
            'Center', 'DisplayName', 'Center', 'Mapping', {'log', ...
            0.002, 0.08}),...
            audioPluginParameter(...
            'pitchChange', 'DisplayName', 'Pitch', 'Mapping', ...
            {'int', -12, 12}));
    end
    
    %----------------------------------------------------------------------
    % PRIVATE PROPERTIES : Used for internal processing and storage 
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % HIDDEN PROPERTIES: Properties accessed outside of plugin, but not
    % visible to user
    %----------------------------------------------------------------------

     properties (Access = private, Hidden)

        % Delay
        pFractionalDelay
       
        % Oscillators
        rModSine
        rModSaw
        rModSquare
        
        
        % mySampleRate Sample rate
        mySampleRate
        
        % Filter
        b = zeros(1,3)
        a = ones(1,3)
        
     end
    
    methods
        
        function p = voiceChangerPlugin

            % LFO for ring Modulation
            
            p.rModSine = audioOscillator('SignalType', 'sine', ...
                    'Frequency', p.finalFreq, 'DCOffset', 0.05, ...
                    'SampleRate', getSampleRate(p));
            p.rModSaw = audioOscillator('SignalType', 'sawtooth', ...
                    'Frequency', p.finalFreq, 'DCOffset', 0.05, ...
                    'Width', 1, 'SampleRate', getSampleRate(p));
            p.rModSquare = audioOscillator('SignalType', 'square', ...
                    'Frequency', p.finalFreq, 'DCOffset', 0.05,...
                    'DutyCycle', 0.5, 'SampleRate', getSampleRate(p));
                
%           dsp.VariableFractionalDelay creates a variable fractional 
%           delay System object that delays a discrete-time input by a 
%           time-varying fractional number of sample periods, as 
%           specified by the second input.
                
            p.pFractionalDelay = dsp.VariableFractionalDelay(...
                'MaximumDelay', 65000);
                
        end
        
       
        
        function out = process(p,in)    
           
                p.finalFreq =  p.rModFreq * (nthroot(2, 12)^p.pitchChange);
                updateFrequency(p, p.finalFreq);
                fs = p.mySampleRate;
                waveTypeDouble = double(p.waveType);
                
                if(length(waveTypeDouble) < 5)
                    oscillator1 = (p.rModSine);
                elseif(length(waveTypeDouble) > 6)
                    oscillator1 =(p.rModSaw);
                elseif(length(waveTypeDouble) < 7)
                    oscillator1 = (p.rModSquare);
                else
                    oscillator1 = p.rModSine;
                end
                                                
                delayInSamples = p.fbDelayNoteDuration * fs;

                
                [d1, d2] = size(in);
                oscillator1.SamplesPerFrame = d1;

                
                delayVector = zeros(d1, d2, 2);
                delayVector(:,:,1) = repmat(delayInSamples+oscillator1(),1,2);
                delayVector(:,:,2) = repmat(delayInSamples+oscillator1(),1,2);
                y = p.pFractionalDelay(in, delayVector);
                oscillator = step(oscillator1);
                
                mix = p.WetDryMix;
                
                out = (in .* ((1-mix)*[0.8 0.8]) + sum(y,3) .* [mix mix]);
                out = out.*(oscillator * [0.2 0.2]);
                out = filter(p.b,p.a,out,p.z);
                out = out .* p.volume / 100.0;
%                 
%                 out = filter(p.b,p.a, (0.2 * oscillator ...
%                     .* (((1-mix) .* in * 0.8) ... 
%                 + (mix.*sum(y,3)))), p.z);
% 
%                 out = out .* p.volume / 100;

        end
        
            function reset(p)
            % Reset sample rate
            Fs = getSampleRate(p);
            p.mySampleRate = Fs;
            
            % Reset oscillators
            p.rModSine.SampleRate = Fs;
            p.rModSaw.SampleRate = Fs;
            p.rModSquare.SampleRate = Fs;
            
            updateFrequency(p, p.finalFreq);
            updateDCOffset(p, p.waveOffset);
            updateControlAmp(p);
            
            % Reset delay
            reset(p.pFractionalDelay); 
            
            p.z = zeros(2);
            [p.b, p.a] = myPassFilters(p.Fc, Fs, p.filtType, p.Q);
            end
            
            
            % THESE ARE LISTENER METHODS APPARENTLY,
            % : Listeners for plugin property/parameter changes
            % (Specified by end-user) Whatever that means...
            function set.Fc(p, Fc)
                
                p.Fc = Fc;
                Fs = getSampleRate(p);
                [p.b, p.a] = myPassFilters(p.Fc, Fs, p.filtType, p.Q);
                
            end
            
            function set.filtType(p, filtType)
                p.filtType = filtType;
                Fs = getSampleRate(p);
                [p.b, p.a] = myPassFilters(p.Fc, Fs, filtType, p.Q);
            end
     
            function set.rModFreq(p, val)
                p.rModFreq = val;
                updateFrequency(p, p.finalFreq) 
            end
            
            function set.waveOffset(p, val)
                p.waveOffset = val;
                updateDCOffset(p, val)
            end
            
            function set.Center(p, val)
                p.Center = val;
                updateDCOffset(p, val)
                updateControlAmp(p);
            end
            
            function set.Range(p, val)
                p.Range = val;
                updateControlAmp(p);
            end
            
    end
    
    methods (Access = private)
        function updateFrequency(p, val)
            p.rModSine.Frequency = val;
            p.rModSaw.Frequency = val;
            p.rModSquare.Frequency = val;
        end
         function updateDCOffset(p,val)
            p.rModSine.DCOffset     	= val;
            p.rModSaw.DCOffset          = val;
            p.rModSquare.DCOffset      	= val;
         end
        
          function updateControlAmp(p)
            controlSignalAmp    = p.Range*p.Center;
            p.rModSine.Amplitude = controlSignalAmp;
            p.rModSaw.Amplitude = controlSignalAmp;
            p.rModSquare.Amplitude = controlSignalAmp;
          end  
          
          function updateFrameSize(p, val)
              p.rModSine.SamplesPerFrame = val;
              p.rModSaw.SamplesPerFrame = val;
              p.rModSquare.SamplesPerFrame = val;
          end
      end
    
        
end

function[b, a] = myPassFilters(Fc, Fs, type, Q)
    w0 = 2*pi*Fc/Fs;
    alpha = sin(w0)/ (2*(2*Q));
    cosw0 = cos(w0);
    norm = 1/(1+alpha);
        
    if(strcmp(type, 'Low'))
        b = (1-cosw0)*norm * [.5 1 .5];
        a = [(1 + alpha) (-2*cosw0) (1 - alpha)];
    elseif(strcmp(type, 'High'))
        b = (1+cosw0)*norm * [.5 -1 0.5];
        a = [(1 + alpha) (-2*cosw0) (1 - alpha)];
    else
    b = 1;
    a = 1;
    end 
end
