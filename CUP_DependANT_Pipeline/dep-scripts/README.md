This repo contains the code sandwiched between the ML model and the client side of the PearlAI website. 
The responsibilities of this repo are:
  - to identify lofs and gofs in patient data
  - to preprocess patient data 
  - to load and run the newest model on patient data
  - to find drugs associated with essential genes found in patient data
  - to find clinical trials associated with those drugs

If you're interested in how it works, take a look at routes.py and patient_processing.py for more information.
