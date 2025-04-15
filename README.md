# Simulation overview
The unconstrained equation of motion for each body of the multibody system are defined within the Component object as the mass, damping, and stiffness matrices, along with the nonlinear vector, generalized forces, and rotational mapping matrix. These make up the general form

$$ M\ddot{q} + D \dot{q} + K q + f^{non} = f^{gen} $$

where $q$ is the vector of generalized coordinates. The Component object can be used for many applications on it's own including simulating the single body, setting initial conditions, plotting states, and checking the conservation of energy.

Once all body equations of motion are defined within their respective Component objects, a ComponentAssembly object must be defined. This gathers the Components and applies the null space method by assembling the full equaitons of motion as defined in the conference paper. The null space transfromation matrix is defined by the user in this object. Once assembled, this object is used to set initial conditions, apply generalized forces, and simulate the entire system with a few short pre-defined method functions.

The following sections go into further detail on the operation of the simulation, including tutorials for key simulation functions.

# Object Oriented Software Design Fundamentals
The simulation is designed with an object oriented programming framework (OOP). This entails combing all commonly used functions and variables that describe a common system into a new user-defined variable type, or class. MATLAB OOP has come a long way in recent years and hase become very powerful. There are a few key OOP concepts that are integral to understanding the operation of this simulation. If you are not familiar with these concepts, I recommend reading the references.

### Classes, Methods, Properties, and Constructors
- https://www.mathworks.com/products/matlab/object-oriented-programming.html
- https://www.mathworks.com/company/technical-articles/introduction-to-object-oriented-programming-in-matlab.html
- https://www.mathworks.com/help/matlab/matlab_oop/create-a-simple-class.html

### Sub-Classes, Super-Classes and Inheritance
- https://www.mathworks.com/help/simscape/lang/subclassing-and-inheritance.html

### Abstract Methods
-https://www.mathworks.com/help/matlab/matlab_oop/abstract-classes-and-interfaces.html

### Defining methods in seperate files
- https://www.mathworks.com/help/matlab/matlab_oop/methods-in-separate-files.html 

# How to Create a New Component
All component objects are subclasses of the Component class, i.e. there will never be a variable in matlab that is of type "Component", just types that are derived from "Component". 

```matlab
classdef MyComponent < Component
...
```

There are several methods that MUST be defined by the user when creating a subclass of "Component". These are shown in the example folder and listed here. These functions can either be listed in the class definition file (MyComponent.m) which is the file containing the required constructor, or they can be defined in seperate files within the class definition folder (@MyComponent) which contains MyComponent.m. These two styles of method function definition can be used within the same class. See https://www.mathworks.com/help/matlab/matlab_oop/methods-in-separate-files.html for more information.
- Compute M
- Compute K
- Compute D
- Compute

# Base Classes

## State

## State Vector

## Component

## Assembly

# How to create a new component

# How to create a new assembly



