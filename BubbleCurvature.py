import numpy as np

# Curvature Calculation
def curvature(x,y):

    dx=np.gradient(x)
    dy=np.gradient(y)

    ddx=np.gradient(dx)
    ddy=np.gradient(dy)

    c1=np.abs(dx*ddy-dy*ddx)
    c2=(dx*dx+dy*dy)**1.5
    cur=c1/c2
    return cur

# Flip in y axis
def FlipYAxAndJoin(x,y):
    tempX=-1*x[1:][::-1]
    tempY=y[1::][::-1]
    x=np.concatenate((tempX,x))
    y=np.concatenate((tempY,y))
    return x,y
