#include "single_pivot.h"
#include "rolling_pivot.h"
//manages the pivoter for the task_rotate_corrs

class Pivot{
    private:
        SinglePivotOfCorrMat* this_pivoter_sp=NULL;
        RollingPivotOfCorrMat* this_pivoter_rp=NULL;
        std::string rotate_type;
    public:
        void setType(std::string input_type){rotate_type=input_type;}
        std::string getType(){return rotate_type;}
        void initiatePivot(TaskHandler& taskhandler, ArgsHandler& xmlin, LogHelper& xmlout, bool& keep_in_task_map){
            if(rotate_type=="SinglePivot"){
                this_pivoter_sp=SinglePivotOfCorrMat::initiateSinglePivot(taskhandler,xmlin,xmlout,keep_in_task_map);
            }
            else if(rotate_type=="RollingPivot"){
                this_pivoter_rp=RollingPivotOfCorrMat::initiateRollingPivot(taskhandler,xmlin,xmlout,keep_in_task_map);
            }
        }
        void checkInitiate(LogHelper& xmlout, XMLHandler& xml_out){
            if( ((this_pivoter_sp==0)&&(rotate_type=="SinglePivot")) || ((this_pivoter_rp==0)&&(rotate_type=="RollingPivot")) ){
               xmlout.output(xml_out);
               throw(std::runtime_error("Could not initiate "+rotate_type));
            }
        }
        void doRotation(uint tmin, uint tmax, char mode, LogHelper& xmllog){
            if(rotate_type=="SinglePivot") this_pivoter_sp->doRotation(tmin,tmax,mode,xmllog);
            else if(rotate_type=="RollingPivot") this_pivoter_rp->doRotation(tmin,tmax,xmllog);
        }
        void writeRotated(uint tmin, uint tmax, const std::string& corrfile, WriteMode overwrite, LogHelper& xmlout, char mode,
                                        char file_format){
            if(rotate_type=="SinglePivot") this_pivoter_sp->writeRotated(tmin,tmax,corrfile,overwrite,xmlout,mode,file_format);
            else if(rotate_type=="RollingPivot") this_pivoter_rp->writeRotated(tmin,tmax,corrfile,overwrite,xmlout);
        }
        void deletePivoter(bool keep){
            if ( (!keep) && (this_pivoter_sp!=0) ) delete this_pivoter_sp; 
            if ( (!keep) && (this_pivoter_rp!=0) ) delete this_pivoter_rp;
        }
        uint getNumberOfLevels(){
            if(rotate_type=="SinglePivot") return this_pivoter_sp->getNumberOfLevels();
            else if(rotate_type=="RollingPivot") return this_pivoter_rp->getNumberOfLevels();
            return 0;
        }
        bool subtractVEV(){
            if(rotate_type=="SinglePivot"){return this_pivoter_sp->subtractVEV();}
            else if(rotate_type=="RollingPivot"){return this_pivoter_rp->subtractVEV();}
            return false;
        }
        GenIrrepOperatorInfo getRotatedOperator(){
            if(rotate_type=="SinglePivot") return this_pivoter_sp->getRotatedOperator();
            else if(rotate_type=="RollingPivot") return this_pivoter_rp->getRotatedOperator();
        }
        void insertEnergyFitInfo(uint level, const MCObsInfo& energyinfo){
            if(rotate_type=="SinglePivot") this_pivoter_sp->insertEnergyFitInfo(level,energyinfo);
            else if(rotate_type=="RollingPivot") this_pivoter_rp->insertEnergyFitInfo(level,energyinfo);
        } 
        void insertAmplitudeFitInfo(uint level, const MCObsInfo& ampinfo){
            if(rotate_type=="SinglePivot") this_pivoter_sp->insertAmplitudeFitInfo(level,ampinfo);
            else if(rotate_type=="RollingPivot") this_pivoter_rp->insertAmplitudeFitInfo(level,ampinfo);
        }
        void reorderLevelsByFitEnergy(LogHelper& xmllog){
            if(rotate_type=="SinglePivot") this_pivoter_sp->reorderLevelsByFitEnergy(xmllog);
            else if(rotate_type=="RollingPivot") this_pivoter_rp->reorderLevelsByFitEnergy(xmllog);
        }
        bool allEnergyFitInfoAvailable(){
            if(rotate_type=="SinglePivot"){return this_pivoter_sp->allEnergyFitInfoAvailable();}
            else if(rotate_type=="RollingPivot"){return this_pivoter_rp->allEnergyFitInfoAvailable();}
            return false;
        }
        bool allAmplitudeFitInfoAvailable(){
            if(rotate_type=="SinglePivot"){return this_pivoter_sp->allAmplitudeFitInfoAvailable();}
            else if(rotate_type=="RollingPivot"){return this_pivoter_rp->allAmplitudeFitInfoAvailable();}
            return false;
        }
        MCObsInfo getEnergyKey(uint level){
            if(rotate_type=="SinglePivot"){return this_pivoter_sp->getEnergyKey(level);}
            else if(rotate_type=="RollingPivot"){return this_pivoter_rp->getEnergyKey(level);}
        }
        MCObsInfo getAmplitudeKey(uint level){
            if(rotate_type=="SinglePivot"){return this_pivoter_sp->getAmplitudeKey(level);}
            else if(rotate_type=="RollingPivot"){return this_pivoter_rp->getAmplitudeKey(level);}
        }
        void computeZMagnitudesSquared(Matrix<MCEstimate>& ZMagSq){
            if(rotate_type=="SinglePivot") this_pivoter_sp->computeZMagnitudesSquared(ZMagSq);
            else if(rotate_type=="RollingPivot") this_pivoter_rp->computeZMagnitudesSquared(ZMagSq);
        }
        const std::set<OperatorInfo>& getOperators(){
            if(rotate_type=="SinglePivot"){return this_pivoter_sp->getOperators();}
            else if(rotate_type=="RollingPivot"){return this_pivoter_rp->getOperators();}
        }
};