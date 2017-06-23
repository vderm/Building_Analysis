# Building_Analysis
Explanation on how to perform and model heat transfer and controls in buildings.

## Introduction
The heat transfer between buildings and their internal and external environments can be modelled as 1-dimensional heat transfer problems. The walls act as barriers to the heat transfer, they resist to it. Although there exist many methods of modelling this heat transfer, a preferred method is to use an electrical engineering analogy: resistors and capacitors. The resistance of walls or windows become resistors; the thermal mass of buildings act as capacitors. Just like in electrical engineering, resistors in series can be combined into an equivalent resistor to simplify the problem.

This repo describes the equations of heat transfer and how I like to set up my model in python. Towards the end, I describe how to calibrate simple resistor-capacitor (RC) models on building data, and finally to use the model for controls using optimization techniques and future predictions (method called: model based predictive controls).

Throughout time, I will be adding other methods and more examples. For example, calibration of models can also be done using a Kalman filter; neural nets can be used as models but then we lose the physical significance of the parameters (black box models). I also need to add in pictures.

Collaboration is warmly welcome.

Vasken Dermardiros
