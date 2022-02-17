import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import onnx
from onnx import version_converter, helper
import onnxmltools
import onnxruntime
from onnxmltools.convert.common.data_types import FloatTensorType, BooleanTensorType, Int32TensorType, DoubleTensorType, Int64TensorType
from onnxmltools.convert.lightgbm.operator_converters.LightGbm import convert_lightgbm  # noqa

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import confusion_matrix, precision_recall_curve, roc_curve, accuracy_score, precision_score, recall_score, f1_score, auc

from skl2onnx import convert_sklearn, update_registered_converter
from skl2onnx.common.shape_calculator import calculate_linear_classifier_output_shapes  # noqa
from skl2onnx.common.data_types import FloatTensorType

import lightgbm as lgb
from lightgbm import LGBMClassifier

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers.experimental import preprocessing
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense

import keras2onnx
import tf2onnx

def plot_metric_curve():
    if MODEL_KEY == 'LightGBM':
        lgb.plot_metric(model)
        plt.savefig('Metric.png')
        plt.close()


def plot_roc_curve(fper, tper):
    plt.plot(fper, tper, color='red', label='ROC')
    plt.plot([0, 1], [0, 1], color='green', linestyle='--')
    plt.xlabel('False Positive Rate = Wrong Matching Rate')
    plt.ylabel('True Positive Rate = Correct Matching Rate')
    plt.title('Receiver Operating Characteristic Curve')
    plt.legend()
    plt.savefig('ROC.png')
    plt.close()

def plot_pr_curve(fper, tper):
    plt.plot(fper, tper, color='red', label='PR')
    plt.xlabel('Recall = Efficiency')
    plt.ylabel('True Positive Rate = Correct Matching Rate')
    plt.title('Precision Recall Curve')
    plt.legend()
    plt.savefig('PR.png')
    plt.close()

def metric_f1(y, y_pred):
    y_hat = np.round(y_pred)
    return 'F1-SCORE', f1_score(y, y_hat), True

def metric_prauc(y, y_pred):
    precision, recall, thresholds = precision_recall_curve(y, y_pred)    
    metric = auc(recall, precision) 
    return 'PR-AUC', metric, True

def get_external_params():

    import yaml
    import argparse

    global MODEL_KEY
    global INPUT_NAME_KEY 
    global LOAD_PARAMS_KEY 
    global CALC_PARAMS_KEY 
    global OBJECT_KEY 
    global HYPER_PARAM_KEY 
    global TRAIN_PARAMS_KEY 
    global SETUP_PARAMS_KEY 
    global OUTPUT_NAME_KEY
    global DEBUG_MODE_KEY

    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    args =  parser.parse_args()    

    with open(args.config) as config:
        obj  = yaml.safe_load(config)

        MODEL_KEY = obj['MODEL']
        INPUT_NAME_KEY = obj['INPUT_NAME']
        LOAD_PARAMS_KEY = obj['LOAD_PARAMS']
        CALC_PARAMS_KEY = obj['CALC_PARAMS']
        OBJECT_KEY = obj['OBJECT']
        HYPER_PARAM_KEY = obj['HYPER_PARAM']
        TRAIN_PARAMS_KEY = obj['TRAIN_PARAMS'] 
        TRAIN_PARAMS_KEY.extend(OBJECT_KEY)
        SETUP_PARAMS_KEY = obj['SETUP_PARAMS']
        OUTPUT_NAME_KEY = obj['OUTPUT_NAME']
        DEBUG_MODE_KEY = obj['DEBUG_MODE']

def preprocessing_data():

    import uproot

    global df
    global df_train
    global df_valid

    global x
    global y
    global x_train
    global y_train
    global x_valid
    global y_valid

    global df_x
    global df_y
    global df_x_train
    global df_y_train
    global df_x_valid
    global df_y_valid

    tree = uproot.open(INPUT_NAME_KEY+':matchTree')

    params = []

    for key in LOAD_PARAMS_KEY:
        row = tree[key].array()
        params.append(row)
    
    array = np.array(params).T

    df = pd.DataFrame(data=array, columns=LOAD_PARAMS_KEY)    
    
    for key in CALC_PARAMS_KEY:    
        if 'Delta' in key: 
            df[key] = df[CALC_PARAMS_KEY[key][0]] - df[CALC_PARAMS_KEY[key][1]]
        elif 'Ratio' in key:
            df[key] = df[CALC_PARAMS_KEY[key][0]] / df[CALC_PARAMS_KEY[key][1]]

    df = df[TRAIN_PARAMS_KEY]
    df_correct = df[df['Truth']==1]
    df_wrong = df[df['Truth']==0]

    df_train, df_valid = train_test_split(df,
                                          test_size=SETUP_PARAMS_KEY['TEST_SIZE'],
                                          random_state=SETUP_PARAMS_KEY['TEST_SEED'])

    df_x_train_correct = df_train[df_train['Truth']==1]
    df_x_train_wrong = df_train[df_train['Truth']==0]
    df_x_train_correct = df_x_train_correct.drop(OBJECT_KEY[0], axis = 1)
    df_x_train_wrong = df_x_train_wrong.drop(OBJECT_KEY[0], axis = 1)
    df_x_train = df_train.drop(OBJECT_KEY[0], axis = 1)
    df_y_train = df_train[OBJECT_KEY[0]]
    
    df_x_valid_correct = df_valid[df_valid['Truth']==1]
    df_x_valid_wrong = df_valid[df_valid['Truth']==0]
    df_x_valid_correct = df_x_valid_correct.drop(OBJECT_KEY[0], axis = 1)
    df_x_valid_wrong = df_x_valid_wrong.drop(OBJECT_KEY[0], axis = 1)    
    df_x_valid = df_valid.drop(OBJECT_KEY[0], axis = 1)
    df_y_valid = df_valid[OBJECT_KEY[0]]

    df_x = df.drop(OBJECT_KEY[0], axis = 1)
    df_y = df[OBJECT_KEY[0]]
        
    x_train = df_x_train.values
    y_train = df_y_train.values

    x_valid = df_x_valid.values
    y_valid = df_y_valid.values

    x = df_x.values
    y = df_y.values

def create_model_TensorFlowNN():
    
    import pickle
    global model

    normalize = preprocessing.Normalization()

    model = tf.keras.Sequential([normalize,
                                 layers.Dense(x_train.shape[1], use_bias=False),
                                 layers.Dense(x_train.shape[1], activation='sigmoid',use_bias=False),
                                 layers.Dense(x_train.shape[1], activation='sigmoid',use_bias=False),
                                 layers.Dense(x_train.shape[1], activation='sigmoid',use_bias=False),
                                 layers.Dense(x_train.shape[1], activation='sigmoid',use_bias=False),
                                 layers.Dense(x_train.shape[1], activation='sigmoid',use_bias=False),
                                 layers.Dense(x_train.shape[1], activation='sigmoid',use_bias=False),
                                 layers.Dense(x_train.shape[1], activation='sigmoid',use_bias=False),
                                 layers.Dense(x_train.shape[1], activation='sigmoid',use_bias=False),
                                 layers.Dense(x_train.shape[1], activation='sigmoid',use_bias=False),
                                 layers.Dense(1, activation='sigmoid')])                             

    model.compile(loss='binary_crossentropy', optimizer=tf.optimizers.Adam(learning_rate=0.005), metrics=['accuracy'])

    training_history = model.fit(x_train, y_train,
                                 epochs=1500,
                                 batch_size=1024,
                                 verbose=1,
                                 validation_data=(x_valid, y_valid))

    with open(OUTPUT_NAME_KEY['ORIGIN'],"wb") as f:
        model.save('tf_model')
        #pickle.dump(model,f)

def convert_onnx_TensorFlowNN():
    
    global onnx_model

    onnx_model, external_tensor_storage = tf2onnx.convert.from_keras(model,input_signature=None, opset=None, custom_ops=None,
                                                                      custom_op_handlers=None, custom_rewriter=None,
                                                                      inputs_as_nchw=None, extra_opset=None, shape_override=None,
                                                                      target=None, large_model=False, output_path=None)
    onnxmltools.utils.save_model(onnx_model,OUTPUT_NAME_KEY['ONNX'])

def create_model_LightGBM():
    
    import pickle

    global model
    
    model = LGBMClassifier(boosting_type=HYPER_PARAM_KEY['TYPE'],
                           objective=HYPER_PARAM_KEY['OBJECTIVE'],
                           learning_rate=HYPER_PARAM_KEY['LEARN_RATE'],
                           max_depth=HYPER_PARAM_KEY['MAX_DEPTH'],
                           n_estimators=HYPER_PARAM_KEY['N_ESTIMATOR'],
                           metric=HYPER_PARAM_KEY['METRIC'])
    
    params = {'boosting_type':HYPER_PARAM_KEY['TYPE'],
              'objective':HYPER_PARAM_KEY['OBJECTIVE'],
              'learning_rate':HYPER_PARAM_KEY['LEARN_RATE'],
              'max_depth':HYPER_PARAM_KEY['MAX_DEPTH'],
              'n_estimators':HYPER_PARAM_KEY['N_ESTIMATOR'],
              'metric':HYPER_PARAM_KEY['METRIC']}
    
    lgb_train = lgb.Dataset(x_train, y_train)
    lgb_valid = lgb.Dataset(x_valid, y_valid, reference=lgb_train)

    model = model.fit(x_train,y_train,
                      eval_metric=[metric_prauc],
                      eval_set=[(x_valid, y_valid),(x_train,y_train)],                      
                      eval_names=['validation','training'])
    
    

    with open(OUTPUT_NAME_KEY['ORIGIN'],"wb") as f:
        pickle.dump(model,f)

def convert_onnx_LightGBM():

    global onnx_model

    update_registered_converter(
        LGBMClassifier, 'LightGbmLGBMClassifier',
        calculate_linear_classifier_output_shapes,convert_lightgbm,
        options={'nocl': [True, False]}
    )

    initial_types = [['inputs', FloatTensorType([None,x_train.shape[1]])]]

    onnx_model = convert_sklearn(model,'lightgbm',initial_types,target_opset=9)

    onnxmltools.utils.save_model(onnx_model,OUTPUT_NAME_KEY['ONNX'])

def get_pred(x):
    
    session_model = onnxruntime.InferenceSession(OUTPUT_NAME_KEY['ONNX'])

    if MODEL_KEY == 'LightGBM':
        y_pred = model.predict(x)
        y_pred_proba = model.predict_proba(x)[:,1]
        input_name = session_model.get_inputs()[0].name
        output_name1 = session_model.get_outputs()[0].name #labels
        output_name2= session_model.get_outputs()[1].name #probabilitoes
        y_pred_onnx = np.squeeze(np.array(session_model.run([output_name1], {input_name: x.astype(np.float32)})))
        y_pred_onnx_proba = np.squeeze(np.array(session_model.run([output_name2], {input_name: x.astype(np.float32)})))

    elif MODEL_KEY == 'TensorFlowNN':
        y_pred = np.round(model.predict(x))
        y_pred_proba = model.predict(x)[:,0]
        input_name = session_model.get_inputs()[0].name
        output_name1 = session_model.get_outputs()[0].name #labels
        y_pred_onnx_proba = session_model.run([output_name1], {input_name: x.astype(np.float32)})[0]
        y_pred_onnx = np.round(y_pred_onnx_proba)

    return y_pred, y_pred_proba, y_pred_onnx, y_pred_onnx_proba

def show_training_valiables(df):
    corr = df.corr()
    sns.heatmap(corr, cmap=sns.color_palette('coolwarm',10))    

def show_results():    
    
    x_check = x
    y_check = y
    df_check = df

    y_pred, y_pred_proba, y_pred_onnx, y_pred_onnx_proba = get_pred(x_check)

    plot_roc_curve(roc_curve(y_check,y_pred_proba)[0],roc_curve(y_check,y_pred_proba)[1])
    plot_pr_curve(precision_recall_curve(y_check,y_pred_proba)[0],precision_recall_curve(y_check,y_pred_proba)[1])
    plot_metric_curve()

    df_check_pred = df_check
    df_check_pred['y_pred'] = np.array(y_pred) 
    
    df_x_check_P = df_check_pred[df_check_pred['Truth']==1]
    df_x_check_N = df_check_pred[df_check_pred['Truth']==0]

    df_x_check_TP = df_check_pred[df_check_pred['y_pred']==1]
    df_x_check_TP = df_x_check_TP[df_x_check_TP['Truth']==1]

    df_x_check_FP = df_check_pred[df_check_pred['y_pred']==1]
    df_x_check_FP = df_x_check_FP[df_x_check_FP['Truth']==0]

    df_x_check_FN = df_check_pred[df_check_pred['y_pred']==0]
    df_x_check_FN = df_x_check_FN[df_x_check_FN['Truth']==1]

    cols = df_x_check_P.columns
    
    axis_scale = 3

    if DEBUG_MODE_KEY['SHOW_PLOT'] == 1:
        for icol in range(len(cols)):
            
            mean_x = df_x_check_N[cols[icol]].mean()
            std_x = df_x_check_N[cols[icol]].std()

            min_range_x = mean_x - axis_scale*std_x
            max_range_x = mean_x + axis_scale*std_x
            
            if cols[icol] == 'Truth' or  cols[icol] == 'y_pred':
                min_range_x = -0.5
                max_range_x =  1.5

            for jcol in range(icol+1,len(cols)):
                
                mean_y = df_x_check_N[cols[jcol]].mean()
                std_y = df_x_check_N[cols[jcol]].std()
                
                min_range_y = mean_y - axis_scale*std_y
                max_range_y = mean_y + axis_scale*std_y

                if cols[jcol] == 'Truth' or  cols[jcol] == 'y_pred':
                    min_range_y = -0.5
                    max_range_y =  1.5

                plt.hist2d(df_x_check_P[cols[icol]], df_x_check_P[cols[jcol]],cmin=1,bins=50,
                           range=[[min_range_x,max_range_x],
                                  [min_range_y,max_range_y]])
                plt.xlabel(cols[icol])
                plt.ylabel(cols[jcol])
                plt.colorbar()
                print('P_Corr_' +cols[icol] + '_' + cols[jcol] + '.png')
                plt.savefig('P_Corr_' +cols[icol] + '_' + cols[jcol] + '.png')
                plt.close()

                plt.hist2d(df_x_check_N[cols[icol]], df_x_check_N[cols[jcol]],cmin=1,bins=50,
                           range=[[min_range_x,max_range_x],
                                  [min_range_y,max_range_y]])
                plt.xlabel(cols[icol])
                plt.ylabel(cols[jcol])
                plt.colorbar()
                print('N_Corr_' +cols[icol] + '_' + cols[jcol] + '.png')
                plt.savefig('N_Corr_' +cols[icol] + '_' + cols[jcol] + '.png')
                plt.close()

                plt.hist2d(df_x_check_TP[cols[icol]], df_x_check_TP[cols[jcol]],cmin=1,bins=50,
                           range=[[min_range_x,max_range_x],
                                  [min_range_y,max_range_y]])
                plt.xlabel(cols[icol])
                plt.ylabel(cols[jcol])
                plt.colorbar()
                print('TP_Corr_' +cols[icol] + '_' + cols[jcol] + '.png')
                plt.savefig('TP_Corr_' +cols[icol] + '_' + cols[jcol] + '.png')
                plt.close()

                plt.hist2d(df_x_check_FP[cols[icol]], df_x_check_FP[cols[jcol]],cmin=1,bins=50,
                           range=[[min_range_x,max_range_x],
                                  [min_range_y,max_range_y]])
                plt.xlabel(cols[icol])
                plt.ylabel(cols[jcol])
                plt.colorbar()
                print('FP_Corr_' +cols[icol] + '_' + cols[jcol] + '.png')
                plt.savefig('FP_Corr_' +cols[icol] + '_' + cols[jcol] + '.png')
                plt.close()

                plt.hist2d(df_x_check_FN[cols[icol]], df_x_check_FN[cols[jcol]],cmin=1,bins=50,
                           range=[[min_range_x,max_range_x],
                                  [min_range_y,max_range_y]])
                plt.xlabel(cols[icol])
                plt.ylabel(cols[jcol])
                plt.colorbar()
                print('FN_Corr_' +cols[icol] + '_' + cols[jcol] + '.png')
                plt.savefig('FN_Corr_' +cols[icol] + '_' + cols[jcol] + '.png')
                plt.close()

    cm = confusion_matrix(y_check, y_pred)
    tn, fp, fn, tp = cm.flatten()
    print('TP  '+str(tp) )
    print('FP  '+str(fp) )
    print('TN  '+str(tn) )
    print('FN  '+str(fn) )
    print(cm)
    
    if tp == 0 and fp == 0:
        print('No positive prediction...')
        exit()

    print("original accuracy : ",accuracy_score(y_check,y_pred))
    print("original precision: ",precision_score(y_check,y_pred))
    print("original recall   : ",recall_score(y_check,y_pred))
    print("original f1       : ",f1_score(y_check,y_pred))
    
    print("onnxruntime accuracy : ",accuracy_score(y_check,y_pred_onnx))
    print("onnxruntime precision: ",precision_score(y_check,y_pred_onnx))
    print("onnxruntime recall   : ",recall_score(y_check,y_pred_onnx))
    print("onnxruntime f1       : ",f1_score(y_check,y_pred_onnx))

    if DEBUG_MODE_KEY['EXPORT_TXT'] == 1:
        export_results_text(y_pred_proba,y_pred_onnx_proba)

def export_results_text(origin, onnx):

    with open('original_probs.txt',"w") as f:
        for p in origin:
            f.write(str(p) + '\n')
        f.close()

    with open('onnx_probs.txt',"w") as f:
        for p in onnx:
            f.write(str(p) + '\n')
        f.close()

def main():
    
    if MODEL_KEY == 'LightGBM': 
        create_model_LightGBM()
        convert_onnx_LightGBM()
    elif MODEL_KEY == 'TensorFlowNN': 
        create_model_TensorFlowNN()
        convert_onnx_TensorFlowNN()

    if DEBUG_MODE_KEY['SHOW_RESULT'] == 1:
        show_results()


if __name__ == "__main__":

    global MODEL_KEY
    global INPUT_NAME_KEY 
    global LOAD_PARAMS_KEY 
    global CALC_PARAMS_KEY 
    global OBJECT_KEY 
    global HYPER_PARAM_KEY 
    global TRAIN_PARAMS_KEY 
    global SETUP_PARAMS_KEY 
    global OUTPUT_NAME_KEY
    global DEBUG_MODE_KEY
    
    global df
    global df_train
    global df_valid

    global x
    global y
    global y_pred
    global y_pred_onnx
    global x_train
    global y_train
    global x_valid
    global y_valid
    global df_x
    global df_y
    global df_x_train
    global df_y_train
    global df_x_valid
    global df_y_valid
    
    global model
    global onnx_model
    
    get_external_params()    
    preprocessing_data()

    main()
