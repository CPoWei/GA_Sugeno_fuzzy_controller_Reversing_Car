# GA reversing car with sugeno(TSK) type fuzzy controller
Train fuzzy controller with Genetic Algorithm for reversing car

Overview : 
1. Here is sugeno type fuzzy control model.
2. The file called “SUGENOwithGA.m” is the main code.
3. Following the main code guide, the GA will start to evolve our fuzzy controller to reach our goal!

Our goal :
1. The goal is reversing the car(represented as a triangle, the sharp angle is the head of the car) to the location around  
   (50,100) with the angle "phi"(calculated from the x-axis to the axis that crosses the head and the tail of the car) around 
   90 degrees.
2. Minimizing "docking_error", which is defined as below: ![alt tag](https://user-images.githubusercontent.com/34533532/34327436-c9bc7c94-e8ff-11e7-98c2-201a4582dbc6.png)

   X_f is equals to 50, Y_f is equals to 100 and Phi_f is equals to 90 degrees.
3. You can find there are a "trajectory_error" in the code, which is only used for seeing how efficiency of the car reversing.![alt tag](https://user-images.githubusercontent.com/34533532/34327437-c9e83424-e8ff-11e7-9faa-c35a8a4bdf0d.png)

Parameter : 
1. X is ranges from 0 to 100
2. Y is ranges from 0 to 100
3. Phi is ranges from -90 to 270 degrees
4. Theta(represented as the angle of the tire can rotate) is ranges from -30 to 30 degrees

Insight : 
1. The probability of mutation should be large. Otherwise, every gene will be same as each other when evolving.
2. Always retaining the outsanding gene seems to be a good way.

Learning Curve : 
![alt tag](https://user-images.githubusercontent.com/34533532/34327437-c9e83424-e8ff-11e7-9faa-c35a8a4bdf0d.png)
