function obj = setDrawingVars(obj)
    %METHOD1 Summary of this method goes here
    %   Detailed explanation goes here
    x_base = [2.3986,2.3986,20.3352,20.3352];
    y_base = [4.7020,0.0500,0.0500,4.7020];

    %corner relative to base vertices
    x_corn = [-2.2494,-1.6933,-1.6672,1.6672,1.6933,2.2494];
    y_corn = [0.0500,0.0500,0.3000,0.3000,0.0500,0.0500];

    %header relative to base verices
    r_head = 2.2500;
    np = 20;
    rad = linspace(0,pi,np);
    x_head = 2.1933*[1,1,-1,-1];
    y_head = [0,0.9055,0.905,0]+0.0500;

    x_head = [x_corn,x_head];
    y_head = [y_corn,y_head];

    obj.x_body = 25.4e-3*[-x_corn+x_base(1),x_corn+x_base(2),x_corn+x_base(3),-x_corn+x_base(4)]';
    obj.y_body = 25.4e-3*[y_corn+y_base(1),-y_corn+y_base(2),-y_corn+y_base(3),y_corn+y_base(4)]';

    obj.x_headers = 25.4e-3*[x_head+x_base(1);...
        -x_head+x_base(3);...
        -x_head+x_base(4);...
        x_head+x_base(2)]';
    obj.y_headers = 25.4e-3*[y_head+y_base(1);...
        -y_head+y_base(3);...
        y_head+y_base(4);...
        -y_head+y_base(2)]';
end
