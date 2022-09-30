function []=MainOptim(o_,o_test,obj_)
    nrun=2;
    OPERATORS={'Alan','Arash','David','DL'};
    GTPARS=[0.6 1.3 1.0 3.6];

    meshmain='/xlwork4/2026068l/PhD/projects/NewParameterization/Inference/GeometryEffect/data/meshes';

    DIRAL=strcat(meshmain,'/Alan');
    DIRAR=strcat(meshmain,'/Arash');
    DIRD=strcat(meshmain,'/David');
    DIRDL=strcat(meshmain,'/DL');
    SUBJECTS=[1009 1015 1042 1079 2008 2024];
    %o_=1; %index of main operator

    MeshesAlan=zeros(6,17376);
    MeshesDavid=zeros(6,17376);
    MeshesArash=zeros(6,17376);
    MeshesDL=zeros(6,17376);

    c=1;
    for subj=SUBJECTS
        load(strcat(DIRAL,'/',num2str(subj),'_1/earlyDiastole/abaqusInputData'))
        mesh=reshape(abaqusInputData.node(:,1:3)',1,17376);
        MeshesAlan(c,:)=mesh;
        load(strcat(DIRAR,'/',num2str(subj),'_1/earlyDiastole/abaqusInputData'))
        mesh=reshape(abaqusInputData.node(:,1:3)',1,17376);
        MeshesArash(c,:)=mesh;
        load(strcat(DIRD,'/',num2str(subj),'_1/earlyDiastole/abaqusInputData'))
        mesh=reshape(abaqusInputData.node(:,1:3)',1,17376);
        MeshesDavid(c,:)=mesh;
        load(strcat(DIRDL,'/',num2str(subj),'_1/earlyDiastole/abaqusInputData'))
        mesh=reshape(abaqusInputData.node(:,1:3)',1,17376);
        MeshesDL(c,:)=mesh;
        c=c+1;
    end


    VolumesAlan=[];
    VolumesArash=[];
    VolumesDavid=[];
    VolumesDL=[];
    cd('/xlwork4/2026068l/PhD/projects/NewParameterization/Inference/GeometryEffect/code')
    for i=1:6
        v=calculate_volumes(MeshesAlan(i,:));
        VolumesAlan=[VolumesAlan v];
        v=calculate_volumes(MeshesArash(i,:));
        VolumesArash=[VolumesArash v];
        v=calculate_volumes(MeshesDavid(i,:));
        VolumesDavid=[VolumesDavid v];
        v=calculate_volumes(MeshesDL(i,:));
        VolumesDL=[VolumesDL v];
    end

    MESHES={MeshesAlan,MeshesArash,MeshesDavid,MeshesDL};
    VOLUMES={VolumesAlan,VolumesArash,VolumesDavid,VolumesDL};
    %%

    %OBJECTIVES={'VS','V','S'};

    MainOperator=OPERATORS{o_};
    MainMeshes=MESHES{o_};
    MainVolumes=VOLUMES{o_};
%     inds=find([1 2 3 4]~=o_);
%     TestMeshes=MESHES(inds);
%     TestVolumes=VOLUMES(inds);

    if ~isfile(strcat('/xlwork4/2026068l/PhD/projects/NewParameterization/Inference/GeometryEffect/data/sims/',MainOperator,'/GTsims.mat'))
        GTsims=zeros(size(MainMeshes,1),25);
        for i=1:size(MainMeshes)
            cd('/xlwork4/2026068l/PhD/Alan_Simulator/LVModelSimulator')
            GTsim=runsim([GTPARS 8],MainMeshes(i,:),strcat('/xlwork4/2026068l/PhD/Alan_Simulator/AbaqusWorkingSpace/',MainOperator),true); %simfun must be defined somewhere
            GTsims(i,:)=GTsim;
        end
        if ~exist(strcat('/xlwork4/2026068l/PhD/projects/NewParameterization/Inference/GeometryEffect/data/sims/',MainOperator),'dir')
            mkdir(strcat('/xlwork4/2026068l/PhD/projects/NewParameterization/Inference/GeometryEffect/data/sims/',MainOperator));
        end
        save(strcat('/xlwork4/2026068l/PhD/projects/NewParameterization/Inference/GeometryEffect/data/sims/',MainOperator,'/GTsims.mat'),'GTsims');
    else
        load(strcat('/xlwork4/2026068l/PhD/projects/NewParameterization/Inference/GeometryEffect/data/sims/',MainOperator,'/GTsims.mat'))
    end
    results.GTpars=GTPARS;
    results.GTsims=GTsims;
    results.MainSimsStandardized=GTsims;
    results.MainSimsStandardized(:,1)=(MainVolumes'-results.MainSimsStandardized(:,1))./results.MainSimsStandardized(:,1);
    %results.MainSimsStandardized(:,1)=results.MainSimsStandardized(:,1)./5;
    results.MainOperator=OPERATORS{o_};
    results.TestOperator=OPERATORS{o_test};
    %%
    resdir=strcat('/xlwork4/2026068l/PhD/projects/NewParameterization/Inference/GeometryEffect/results/',obj_,num2str(nrun));
    if ~exist(resdir,'dir')
            mkdir(resdir);
    end

    results.optims=zeros(6,4);
    for i_=1:numel(SUBJECTS) %Each individual subject
        tracedir=strcat(resdir,'/lossTrace',results.MainOperator,'-',results.TestOperator,num2str(i_),'.txt');
%        for operator=1:3 %different operators
        CurrentTestMeshes=MESHES{o_test};
        CurrentTestVolumes=VOLUMES{o_test};
        GTsim=results.MainSimsStandardized(i_,:);
        %for obj_=1:3 %V,S or VS
        cd('/xlwork4/2026068l/PhD/projects/NewParameterization/Inference/GeometryEffect/code');

        obj=@(x) loss(x,GTsim,CurrentTestMeshes(i_,:),...
            CurrentTestVolumes(i_),results.MainOperator,results.TestOperator,...
            obj_);
        options= optimoptions('fmincon');
        options.Algorithm='sqp';
        options.DiffMinChange=1e-4;
        options.DiffMaxChange=1e-3;
        options.OutputFcn=@(x,y,z) outputfn(x,y,z,tracedir);
        options.OptimalityTolerance=1e-6;
        x = fmincon(obj,GTPARS,[],[],[],[],[0.05 0.05 0.05 0.05],[10 10 10 10],[],options);
        results.optims(i_,:)=x;
        cd('/xlwork4/2026068l/PhD/projects/NewParameterization/Inference/GeometryEffect/code');
        
        save(strcat(resdir,'/',results.MainOperator,'-',results.TestOperator),'results')
        rmdir(strcat('/xlwork4/2026068l/PhD/Alan_Simulator/AbaqusWorkingSpace/OPTIM',results.MainOperator,'-',results.TestOperator,obj_),'s')
    end
end
