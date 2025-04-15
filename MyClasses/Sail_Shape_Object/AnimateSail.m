function AnimateSail(obj, t, tip, z_scaling, folder, fileName)
% % Animates the four beam solar sail without sail material. Currently
% % replicating the single current beam 4 times.
% %
% % INPUTS
% % t : 1 x n time vector
% % u : n x 4 array of tip deflections
% % z_scaling : scaling to apply to the z axis of the plottings
% %
% %-------------------------------------------------------------------------
% 
% % Initialize video to save animation
%             % filename = [folder, '/', fileName];
%             % myVideo = VideoWriter(filename, 'MPEG-4'); %open video file
%             % myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
%             % open(myVideo)
% 
% % change in tip deflection
% tip_0 = tip(1,:);
% tip_delta = tip(2:end, :) - tip(1:end-1, :);
% 
% % create the new figue
% myVideoFigHandle = figure();
% 
% % Loop through all time steps
% for lv1 = 1:length(t)
% 
%     % Get current theta and deflection
%     if lv1 ==1 % Nominal tip deflection
%         tip_i = tip_0;
%     else % Update tip deflection to a new value with the delta
%         tip_i = tip_delta(lv1-1,:);
%     end
% 
%     % Pause to let the figure catch up
%     pause(1/100);
% 
%     % Apply the new tip deflection 
%     obj = obj.applyTipDeflection(tip_i(1), tip_i(2), tip_i(3), tip_i(4));
% 
%     %
%     obj.plotSailMeshOnFigure(obj, z_scaling, myVideoFigHandle)
%     ax = gca;
%     ax.XAxis.Visible = 'off';
%     ax.YAxis.Visible = 'off';
%     ax.ZAxis.Visible = 'off';
%     % gcf
% 
%     pause(1/100); %Pause and grab frame
%             % frame = getframe(myVideoFigHandle); %get frame
%             % writeVideo(myVideo, frame);
% 
% end
%             % close(myVideo)
end








