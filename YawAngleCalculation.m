
clc;
clear;
close all; 
%% Step0 随机数构建数据集
        ACC_X_1 = randperm(2001,2000)-1000;        % X轴向加速度
        ACC_Y_1 = (randperm(2002,2000)-1000).*0.2; % Y轴向加速度
        ACC_Z_1 = randperm(2003,2000)-1000;        % Z轴向加速度 
        Rot_X_1 = randperm(2004,2000).*0.05;        % X轴转速
        Rot_Y_1 = randperm(2005,2000).*0.05;        % Y轴转速
        Rot_Z_1 = randperm(2006,2000).*0.05;        % Z轴转速     
    
     % 此时人为设定车平面分布于ZX平面，此情况仅供测试使用
         % 画图
         figure('Name','Fig1_Random_Points','NumberTitle','off');
         scatter3(ACC_X_1,ACC_Y_1,ACC_Z_1);
         view(40,35)
         xlim([-1000 1000]); 
         ylim([-1000 1000]);
         zlim([-1000 1000]);
               xlabel('x axis','FontName','Times newman','FontSize',12);
               ylabel('y axis','FontName','Times newman','FontSize',12);
               zlabel('z axis','FontName','Times newman','FontSize',12);
%% Step1 数据采集   
     % 实际情况下此部分应为实时采集的加速度和转速数据组
     % Time = linspace(1, 2000 ,2000); %此处为手动添加数据点采集对应的时间     
     Time = linspace(1,length(ACC_X_1) ,length(ACC_X_1)); %此处为手动添加数据点采集对应的时间
     Data = [Time',ACC_X_1',ACC_Y_1',ACC_Z_1',Rot_X_1',Rot_Y_1',Rot_Z_1'];  %数据组 采集时间+X轴加速度+Y轴加速度+Z轴加速度
     
%% Step2 数据预处理   
     
     % % 轴向速度判定控制分段，此处省略步骤，实际使用时具体分析情况     
       
%% Step3 利用最小二乘法获得车辆的真实水平面 
 % 3.1 拟合前数据处理
     mean_ACC_all = mean(Data(:,2:4));   % XYZ轴向加速度平均值 mean
     std_ACC_all = std(Data(:,2:4));     % XYZ轴向加速度标准差 σ
     
    % Beta_i 对应了数据组的 Z_score 标准化分值    β_i=(x_i-mean)/σ
     Beta_all = (Data(:,2:4)-mean_ACC_all)./(std_ACC_all);  Beta_all = [Data(:,1),Beta_all];
     % Beta_all_max = max(Beta_all);
         FAR = 5; % FalseAcceptRate 错误接受率 可根据实际情况修改

    Data_Exclude_Outlier = Data (abs(Beta_all(:,2)) < FAR,:);                 % X轴晒选
    Data_Exclude_Outlier = Data_Exclude_Outlier(abs(Beta_all(:,3)) < FAR,:);  % Y轴晒选
    Data_Exclude_Outlier = Data_Exclude_Outlier(abs(Beta_all(:,4)) < FAR,:);  % Z轴晒选

 % 3.2 根据最小二乘法拟合车平面    （此处实际使用了奇异值分解（SVD）的手段求最小二乘平面）
    
         % 三维空间中的二维平面可以表示为 Ax+By+Cz+D = 0
         figure('Name','Fig2_Random_Points_Fit','NumberTitle','off');
         scatter3(ACC_X_1,ACC_Y_1,ACC_Z_1,'r');% scatter3(ACC_X_1,ACC_Y_1,ACC_Z_1,'r','filled');
         view(40,35)
         xlim([-1000 1000]); 
         ylim([-1000 1000]);
         zlim([-1000 1000]);
               xlabel('x axis','FontName','Times newman','FontSize',12);
               ylabel('y axis','FontName','Times newman','FontSize',12);
               zlabel('z axis','FontName','Times newman','FontSize',12);
         
           hold on;

        planeData= Data_Exclude_Outlier(:,2:4);   % 导入参与拟合平面的xyz三列数据

        % 协方差矩阵的SVD变换中，最小奇异值对应的奇异向量就是平面的方向
        xyz0 = mean(planeData,1);                         % 计算列均值
        centeredPlane = bsxfun(@minus,planeData,xyz0);    % bsxfun(@minus,planeData,xyz0)指planeData矩阵中元素与列均值xyz0的差值
        Plane4SVD = centeredPlane' * planeData;  
        [SVD_U,SVD_S,SVD_V]=svd(Plane4SVD);               % SVD分解求特征值 U·S·V

        % 矩阵计算结果中提取二维平面系数
        Plane.A=SVD_V(1,3);         % Ax+By+Cz+D = 0 中系数A
        Plane.B=SVD_V(2,3);         % Ax+By+Cz+D = 0 中系数B
        Plane.C=SVD_V(3,3);         % Ax+By+Cz+D = 0 中系数C
        Plane.D=-dot([Plane.A Plane.B Plane.C],xyz0);  % Ax+By+Cz+D = 0 中系数D 通过ABC构成的向量与列均值构成向量的点乘得到

        % 画图验证拟合结果 实际应用中可不使用 
            fit.x = -1000:20:1000;
            fit.y = -1000:20:1000;
            [FIT.X,FIT.Y]= meshgrid (fit.x,fit.y);
            FIT.Z = -(Plane.D + Plane.A * FIT.X + Plane.B * FIT.Y)/Plane.C;  % Ax+By+Cz+D = 0
            mesh(FIT.X,FIT.Y,FIT.Z);
            hold on 
            FIT.Z2 = zeros(size(FIT.Z));
            mesh(FIT.X,FIT.Y,FIT.Z2);
            
     
%% Step4：根据车平面对实际采集数据组进行轴角转换    
      Normal_Vector_1 = [0,0,1];                           % 初始三维平面法向量  N
      Normal_Vector_CarPlane = [Plane.A,Plane.B,Plane.C];  % 现在的车平面法向量  M
      % 计算初始平面与车平面间夹角 （通过向量的点乘实现）
      cos_RotateTheta =  dot(Normal_Vector_1,Normal_Vector_CarPlane) ./ ( norm(Normal_Vector_1)* norm(Normal_Vector_CarPlane));
      % 计算坐标转换所需的旋转轴向量 （通过向量的叉乘实现）
      Vector_RotateAxis = (cross(Normal_Vector_1,Normal_Vector_CarPlane))/(norm(cross(Normal_Vector_1,Normal_Vector_CarPlane)));
     
      % 计算转动矩阵 Rotate_Matrix （Q）
      c = cos_RotateTheta; s = sqrt(1-c^2); C = 1-c;
        x = Vector_RotateAxis(1,1) ; y = Vector_RotateAxis(1,2) ; z = Vector_RotateAxis(1,3);
       % z = Plane.A ; y = Plane.B ; x = Plane.C;
        Q = [x*x*C+c,x*y*C-z*s,x*z*C+y*s;y*x*C+z*s,y*y*C+c,y*z*C-x*s;z*x*C-y*s,z*y*C+x*s,z*z*C+c] ;
    % 对轴向加速度数据进行旋转操作  
        Data_Acc = Data(:,2:4); % Data3 = Data(:,2:4)';

                         % Data_Rotate = dot(Q,(Data(:,2:4)'));  错误
            %  Data_Rotate = Q*Data3;
        Data_Acc_Rotated = Data_Acc*Q;
       
         figure('Name','Fig3_Points_Rotated','NumberTitle','off');
          Data_Rotate2= Data_Acc_Rotated;
        %  scatter3(Data_Rotate(1,:),Data_Rotate(2,:),Data_Rotate(3,:));    
          scatter3(Data_Rotate2(:,1),Data_Rotate2(:,2),Data_Rotate2(:,3));    
         view(40,35)
         xlim([-1000 1000]); 
         ylim([-1000 1000]);
         zlim([-1000 1000]);
               xlabel('x axis','FontName','Times newman','FontSize',12);
               ylabel('y axis','FontName','Times newman','FontSize',12);
               zlabel('z axis','FontName','Times newman','FontSize',12);
    % 对转速数据进行旋转操作           
        Data_Rot = Data(:,5:7);         
        Data_Rot_Rotated = Data_Rot*Q; 
        
        Data_All_Rotated = [Data(:,1),Data_Acc_Rotated,Data_Rot_Rotated];  % 信号采集时间 三轴加速度 三轴转速
         
 %% Step5：根据新数据组和实测的偏航角对照，使用机器学习获得关联性        
       % XY平面上的  1.平均加速度，2.平均转速，3.加速度的标准差，4.转速的标准差 共8个数据构建1×8的矢量
      %  i 时刻
       Data_All_Rotated_i_1 = mean(Data_All_Rotated);  
       Data_All_Rotated_i_2 = std(Data_All_Rotated);
       Data_All_Rotated_i = [Data_All_Rotated_i_1(1,2:3),Data_All_Rotated_i_2(1,2:3),0,0,Data_All_Rotated_i_1(1,5:6),Data_All_Rotated_i_2(1,5:6),0,0];
%         accMeanX 0.10108176
%         accMeanY 0.17573071
%         accStdX 0.00604311
%         accStdY 0.0041152
%         accErrX 0.0             %此处系数为文章中提供，实际需使用机器学习通过实际数据集进行拟合迭代优化结果
%         accErrY 0.0
%         rotMeanX 0.11514755
%         rotMeanY 0.58944343
%         rotStdX 0.00560081
%         rotStdY 0.00283743
%         rotErrX 0.0
%         rotErrY 0.0
       
       feature_importances =[0.10108176;0.17573071;0.00604311;0.0041152;0.0;0.0;0.11514755;0.58944343;0.00560081;0.00283743;0.0;0.0];
       
       % 航向角
       Yaw_Angle_i = Data_All_Rotated_i *  feature_importances;
         
         

         
         
         
         
      
      
     