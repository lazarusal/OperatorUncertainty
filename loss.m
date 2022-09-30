function [fx]=loss(x,y,G,v0,MainOperator,TestOperator,objtype)
    cd('/xlwork4/2026068l/PhD/Alan_Simulator/LVModelSimulator')
    sim=runsim([x 8],G,strcat('/xlwork4/2026068l/PhD/Alan_Simulator/AbaqusWorkingSpace/OPTIM',MainOperator,'-',TestOperator,objtype),true); %simfun must be defined somewhere
    sim(1)=(v0-sim(1))/sim(1);
    if strcmp(objtype,'V')
        fx=(sim(1)-y(1))^2;
    elseif strcmp(objtype,'S')
        fx=sum((sim(2:end)-y(2:end)).^2);
    else
        fx=sum((sim-y).^2);
    end        
    cd('/xlwork4/2026068l/PhD/projects/NewParameterization/Inference/GeometryEffect/code')
end
    