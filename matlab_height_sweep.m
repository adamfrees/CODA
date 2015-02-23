mphstart(2036)
import com.comsol.model.*
import com.comsol.model.util.*
model = mphload('/home/frees/Dropbox/_UW/Scale-up/Comsol_test_1/var_height.mph');

order = [3,1,5,14,15,6,7,8,9,19,10,11,16,12,17];
sweepnum = 1;
looping = true;

target = [1,1,0,0,0.01,0.01,0.01]
numHeights = 11
heightList = linspace(0.08,0.09,numHeights)

%Cheap hack:
%startingVs = [0.3988,0.3482,0.3947,-0.3462,0.1335,-0.3943,-0.1078,-0.01,-0.3577,0.122,-0.3341,-0.101,-0.0841,-0.1032,-0.1]
[dim, loopLength] = size(order);
for i=1:loopLength
    n=order(i);
    model.physics('es').feature(strcat('pot',int2str(n))).set('V0',startingVs(i));
end


for iter=1:numHeights
    count = 0;
    height = heightList(iter);
    model.geom('geom1').feature('wp4').set('quickz',height)
    while true
        try
            model.mesh('mesh1').run()
            break
        catch e
            height = height+0.00001
            model.geom('geom1').feature('wp4').set('quickz',height)
        end
    end
    while looping
        [dim, loopLength] = size(order);
        voltList = [];
        for i =1:loopLength
            n=order(i);
            voltList = [voltList,str2double(model.physics('es').feature(strcat('pot',int2str(n))).getString('V0'))];
        end
        
        
        for i = 0:loopLength
            n = 0;
            if i ~= 0;  n=order(i); end
            if i~=0; volts = str2double(model.physics('es').feature(strcat('pot',int2str(n))).getString('V0')); end
            if i~=0; model.physics('es').feature(strcat('pot',int2str(n))).set('V0',volts-0.0005); end
            if i~=0; voltList(i) = volts-0.005; end
            fprintf(fopen(strcat('~/Dropbox/_UW/Scale-up/Var_height/',int2str(iter),'sweep',int2str(sweepnum),'.txt'),'a+'),'%g ',voltList);
            if i~=0; voltList(i) = volts; end
            model.study('std1').run()
            if i~=0; model.physics('es').feature(strcat('pot',int2str(n))).set('V0',volts); end
            
            model.result.export('data1').set('filename',strcat('/home/frees/Dropbox/_UW/Scale-up/Var_height/',int2str(iter),'run',int2str(count),'.txt'));
            model.result.export('data1').run;
            [out,cmdout] = system(strcat('python /home/frees/Dropbox/_UW/Scale-up/NewSweep/ReadOutput.py ',strcat(' /home/frees/Dropbox/_UW/Scale-up/Var_height/',int2str(iter),'run',int2str(count),'.txt')))
            idx = find(ismember(cmdout,')'),1,'last');
            if cmdout(idx)==')'; cmdout(1:idx)=[]; end
            fprintf(fopen(strcat('~/Dropbox/_UW/Scale-up/Var_height/',int2str(iter),'sweep',int2str(sweepnum),'.txt'),'a+'),'%g ',str2num(cmdout));
            fprintf(fopen(strcat('~/Dropbox/_UW/Scale-up/Var_height/',int2str(iter),'sweep',int2str(sweepnum),'.txt'),'a+'),'%s\n','');
            count = count+1;
            if i==0
                str2num(cmdout)
                nextStep = target-str2num(cmdout);
                nextStep(3)=0;
                nextStep(4)=0;
                nextStep = num2str(nextStep);
                while isempty(strfind(nextStep,'  '))==0; nextStep = strrep(nextStep,'  ',' '); end
                nextStep = strrep(nextStep,' ',',')
            end
        end
        [out,cmdout] = system(strcat(strcat('python /home/frees/Dropbox/_UW/Scale-up/Var_height/findNewPoint.py ',strcat(' /home/frees/Dropbox/_UW/Scale-up/Var_height/',int2str(iter),'sweep',int2str(sweepnum),'.txt')),[' ',nextStep]))
        idx = find(ismember(cmdout,')'),1,'last');
        if cmdout(idx)==')'; cmdout(1:idx)=[]; end
        myStart=2;
        voltChange = []
        idx = find(ismember(cmdout,','),loopLength,'last');
        for i = 1:loopLength-1
            voltChange = [voltChange,str2double(cmdout(myStart:idx(i)-1))];
            myStart = idx(i)+1;
        end
        voltChange = [voltChange,str2double(cmdout(myStart:end))];
        voltList = voltList + voltChange
        for i =1:loopLength
            n=order(i);
            volts = voltList(i)
            model.physics('es').feature(strcat('pot',int2str(n))).set('V0',volts);
        end
        if sweepnum==8;
            looping=false;
        end
        sweepnum = sweepnum+1;
    end
end
