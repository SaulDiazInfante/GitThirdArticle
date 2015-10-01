cat ./Preamble.tex > StrongConvergenceLSMethod.tex
echo '\n%********************************************************************************************' >> StrongConvergenceLSMethod.tex 
echo '%              Section 1' >> StrongConvergenceLSMethod.tex
echo '%********************************************************************************************\n' >> StrongConvergenceLSMethod.tex
echo '\section{Introduction} \n' >> StrongConvergenceLSMethod.tex
cat  ./Introduction.tex  >> StrongConvergenceLSMethod.tex
echo '\n%********************************************************************************************' >> StrongConvergenceLSMethod.tex 
echo '%              Sectison 2' >> StrongConvergenceLSMethod.tex
echo '%********************************************************************************************\n' >> StrongConvergenceLSMethod.tex
echo '\section{General Settings} \n' >> StrongConvergenceLSMethod.tex
cat  ./GeneralSettings.tex >> StrongConvergenceLSMethod.tex
echo '\n%********************************************************************************************' >> StrongConvergenceLSMethod.tex 
echo '%              Section 3' >> StrongConvergenceLSMethod.tex
echo '%********************************************************************************************\n' >> StrongConvergenceLSMethod.tex
echo '	\section{Construction of the Linear Steklov method}  \n' >> StrongConvergenceLSMethod.tex
cat  ./ConstructionoftheLinearSteklovmethod.tex >> StrongConvergenceLSMethod.tex
echo '\n%********************************************************************************************' >> StrongConvergenceLSMethod.tex 
echo '%              Section 4' >> StrongConvergenceLSMethod.tex
echo '%********************************************************************************************\n' >> StrongConvergenceLSMethod.tex
echo '\section{Strong Convergence of the Linear Steklov Method} \n' >> StrongConvergenceLSMethod.tex
cat  ./StrongConvergenceoftheLinearSteklovMethod.tex  >> StrongConvergenceLSMethod.tex
echo '\n%********************************************************************************************' >> StrongConvergenceLSMethod.tex 
echo '%              Section 5' >> StrongConvergenceLSMethod.tex
echo '%********************************************************************************************\n' >> StrongConvergenceLSMethod.tex
echo '	\section{Convergence Rate}' >> StrongConvergenceLSMethod.tex
cat  ./ConvergenceRate.tex>> StrongConvergenceLSMethod.tex
echo '\n%********************************************************************************************' >> StrongConvergenceLSMethod.tex 
echo '%              Section 6' >> StrongConvergenceLSMethod.tex
echo '%********************************************************************************************\n' >> StrongConvergenceLSMethod.tex
echo '	\section{Numerical Simulation}\n' >> StrongConvergenceLSMethod.tex
cat ./NumericalSimulation.tex>> StrongConvergenceLSMethod.tex
echo '\n%********************************************************************************************' >> StrongConvergenceLSMethod.tex 
echo '%\t\t\t\t Section 7' >> StrongConvergenceLSMethod.tex
echo '%********************************************************************************************\n' >> StrongConvergenceLSMethod.tex
echo '\section{Conclusions} \n' >> StrongConvergenceLSMethod.tex
cat ././Conclusions.tex >> StrongConvergenceLSMethod.tex
echo '\n%********************************************************************************************' >> StrongConvergenceLSMethod.tex 
echo '%\t\t\t\t Bibliography' >> StrongConvergenceLSMethod.tex
echo '%********************************************************************************************\n' >> StrongConvergenceLSMethod.tex
#
cat ./BackMater.tex >> StrongConvergenceLSMethod.tex
cat ./SecondArticle.bbl >> StrongConvergenceLSMethod.tex 
echo '\n%********************************************************************************************' >> StrongConvergenceLSMethod.tex 
echo '%\t\t\t\t Apendix' >> StrongConvergenceLSMethod.tex
echo '%********************************************************************************************\n' >> StrongConvergenceLSMethod.tex
echo '\end{document}' >> StrongConvergenceLSMethod.tex
