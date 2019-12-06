# SUBTRwMNPO
Frequency Dependent Subtraction Method with Minimum Norm Projection Operator for Removal of Interference from Fetal MEG Data

Background: A frequency dependent subtraction method, SUBTR, is developed to remove maternal and fetal magnetocardiography (mMCG and fMCG) interference from fetal magnetoencephalography (fMEG). But channels close to fetal head cannot be used as references for SUBTR in order to protect fMEG from subtraction and this results in cardiac residual when these channels have important fMCG frequency components. Cardiac residual creates noise in Evoked Response (ER) which results in poor ER detection. 

New Method:  We developed an enhanced SUBTR algorithm, which we call SUBTR with Minimum Norm Projection Operator (SUBTRwMNPO), by employing covariance based Minimum Norm Projection Operators (MNPO). mMCG and fMCG signals are extracted from the raw data using MNPO and they are subtracted in the frequency domain from raw data to extract fetal Evoked Response (fER).

Results:  When tested on 87 datasets, SUBTRwMNPO is shown to attenuate cardiac interference almost totally resulting in a clean fER signal.

Comparison with Existing Methods: Cardiac attenuation with SUBTRwMNPO is either as good as or better than SUBTR. SUBTRwMNPO has higher attenuation rate for the datasets where SUBTR leaves cardiac residual.

Conclusions: SUBTRwMNPO is successful in removing cardiac interference regardless of the orientation of fMCG and fMEG signal spaces. It can also be used to remove cardiac interference when there is no prior knowledge of fetal head location.

For detailed information please check the paper: " Bisgin, Neslihan; Wilson, James D.; Eswaran, Hari; An Enhanced Frequency Dependent Subtraction Algorithm for Removal of Cardiac Interference from fMEG by Using Minimum Norm Projection Operator
Journal of Neuroscience Methods, 2019"
