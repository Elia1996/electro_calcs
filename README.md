# electro_calcs
## cross_talk.py script
The main implementation in this script is the following class:

### xtalk_circuit
      Returned list:
            [Mne_ind   Mne_cap  Mfe_ind  Mfe_cap  dominating_effect]

        This function take as reference the following circuit:

                                Generator Conductor
            ┌----/\/\/\------------------------>--------------------------┐
            |      Rs                         Ig(z,f)       +             |
            |                     Receptor Conductor      Vg(z,f)         |
            |          ┌------------>----------------------------┐        |
    Vs(f) / + \        |  +      Ir(z,f)   +                  +  |        |
          \_-_/        \                                         /        \
            |      Rne /  Vne            Vr(z,f)             Vfe \ Rfe    / Rl
            |          \                                         /        \
            |          |  -                -      Ig+Ir       -  |        |
            └----------┴---------------------------->------------┴--------┘
                                Refecence conductor
                        |-----------------------------------------|
                                            L
                        |-----------------------------------------|---------------------> z
                        z=0                                       z=L

