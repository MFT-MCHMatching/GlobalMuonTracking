import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

import onnx
from onnx import version_converter, helper
import onnxmltools
import onnxruntime
from onnxmltools.convert.common.data_types import FloatTensorType, BooleanTensorType, Int32TensorType, DoubleTensorType, Int64TensorType
from onnxmltools.convert.lightgbm.operator_converters.LightGbm import convert_lightgbm  # noqa

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score

from skl2onnx import convert_sklearn, update_registered_converter
from skl2onnx.common.shape_calculator import calculate_linear_classifier_output_shapes  # noqa
from skl2onnx.common.data_types import FloatTensorType

import lightgbm as lgb
from lightgbm import LGBMClassifier

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
#from tensorflow.keras.layers.experimental import preprocessing
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense

import keras2onnx
import tf2onnx

#from sklearn.metrics import auc
#from sklearn.metrics import precision_recall_curve

def getExternalParams():

    import yaml
    import argparse

    global model_key
    global infile_key 
    global load_key 
    global calc_key 
    global object_key 
    global hyper_param_key 
    global df_key 
    global train_params_key 
    global output_name_key
    global debug_mode_key

    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    args =  parser.parse_args()    

    with open(args.config) as config:
        obj  = yaml.safe_load(config)

        model_key = obj['MODEL']
        infile_key = obj['INPUT_NAME']
        load_key = obj['LOAD_PARAMS']
        calc_key = obj['CALC_PARAMS']
        object_key = obj['OBJECT']
        hyper_param_key = obj['HYPER_PARAM_LIGHTGBM']
        df_key = obj['TRAIN_PARAMS'] 
        df_key.extend(object_key)
        train_params_key = obj['SETUP_PARAMS']
        output_name_key = obj['OUTPUT_NAME']
        debug_mode_key = obj['DEBUG_MODE']

def preprocessingData():

    import uproot

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

    tree = uproot.open(infile_key+':matchTree')

    params = []

    for key in load_key:
        row = tree[key].array()
        params.append(row)
    
    array = np.array(params).T

    df = pd.DataFrame(data=array, columns=load_key)    
    
    for key in calc_key:    
        if 'Delta' in key: 
            df[key] = df[calc_key[key][0]] - df[calc_key[key][1]]
        elif 'Ratio' in key:
            df[key] = df[calc_key[key][0]] / df[calc_key[key][1]]

    df = df[df_key]

    df_train, df_valid = train_test_split(df,
                                          test_size=train_params_key['TEST_SIZE'],
                                          random_state=train_params_key['TEST_SEED'])

    df_x_train = df_train.drop(object_key[0], axis = 1)
    df_y_train = df_train[object_key[0]]

    df_x_valid = df_valid.drop(object_key[0], axis = 1)
    df_y_valid = df_valid[object_key[0]]

    df_x = df.drop(object_key[0], axis = 1)
    df_y = df[object_key[0]]
        
    x_train = df_x_train.values
    y_train = df_y_train.values

    x_valid = df_x_valid.values
    y_valid = df_y_valid.values

    x = df_x.values
    y = df_y.values

def f1_metric(preds, train_data):
    labels = train_data.get_label()
    return 'f1', f1_score(labels, preds, average='weighted'), True

def train_TensorFlowNN():

    import pickle
    global model

    model = tf.keras.Sequential([#normalize,
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

    model.compile(loss='binary_crossentropy', optimizer=tf.optimizers.Adam(learning_rate=0.0005), metrics=['accuracy'])

    model.summary

    training_history = model.fit(x_train, y_train,
                                 epochs=150,
                                 batch_size=1000,
                                 verbose=1,
                                 validation_data=(x_eval, y_eval))

    with open(output_name_key['ORIGIN'],"wb") as f:
        pickle.dump(model,f)

def convert_ONNX_TensorFlowNN():

    global onnx_model

    onnx_model, external_tensor_storage = tf2onnx.convert.from_keras(model,input_signature=None, opset=None, custom_ops=None,
                                                                      custom_op_handlers=None, custom_rewriter=None,
                                                                      inputs_as_nchw=None, extra_opset=None, shape_override=None,
                                                                      target=None, large_model=False, output_path=None)
    onnxmltools.utils.save_model(onnx_model,output_name_key['ONNX'])


def train_LightGBM():
    
    import pickle
    global model

    model = LGBMClassifier(boosting_type=hyper_param_key['TYPE'],
                           objective=hyper_param_key['OBJECTIVE'],
                           learning_rate=hyper_param_key['LEARN_RATE'],
                           max_depth=hyper_param_key['MAX_DEPTH'],
                           n_estimators=hyper_param_key['N_ESTIMATOR'],
                           metric=hyper_param_key['METRIC'])
    
    training_history = model.fit(x_train,y_train,eval_metric=hyper_param_key['METRIC'],
                                 eval_set=[(x_train, y_train),
                                           (x_valid, y_valid)],
                                 eval_names=['train', 'validation'])
    
    lgb.plot_metric(model)

    with open(output_name_key['ORIGIN'],"wb") as f:
        pickle.dump(model,f)

def convert_ONNX_LightGBM():

    global onnx_model

    update_registered_converter(
        LGBMClassifier, 'LightGbmLGBMClassifier',
        calculate_linear_classifier_output_shapes,convert_lightgbm,
        options={'nocl': [True, False]}
    )

    initial_types = [['inputs', FloatTensorType([None,x_train.shape[1]])]]

    onnx_model = convert_sklearn(model,'lightgbm',initial_types,target_opset=9)

    onnxmltools.utils.save_model(onnx_model,output_name_key['ONNX'])

def showFeature(df):
    corr = df.corr()
    sb.heatmap(corr, cmap=sb.color_palette('coolwarm',10))    

def showAccuracyPrecision():
    pred_model = model.predict(x)
    pred_model_proba = model.predict_proba(x)
    session_model = onnxruntime.InferenceSession(output_name_key['ONNX'])

    input_name = session_model.get_inputs()[0].name
    output_name1 = session_model.get_outputs()[0].name #labels
    output_name2= session_model.get_outputs()[1].name #probabilitoes

    pred_onnx_model = np.squeeze(np.array(session_model.run([output_name1], {input_name: x.astype(np.float32)})))
    pred_onnx_model_proba = np.squeeze(np.array(session_model.run([output_name2], {input_name: x.astype(np.float32)})))

    print("original accuracy:  ",accuracy_score(y,pred_model))
    print("original precision: ",precision_score(y,pred_model))

    print("onnxruntime accuracy:  ",accuracy_score(y,pred_onnx_model))
    print("onnxruntime precision: ",precision_score(y,pred_onnx_model))

    if debug_mode_key['EXPORT_TXT'] == 1:
        exportProbsAsText(pred_model_proba,pred_onnx_model_proba)

def exportProbsAsText(origin, onnx):

    with open('original_probs.txt',"w") as f:
        for p in origin:
            _p = float(p[1])
            f.write(str(_p) + '\n')
        f.close()

    with open('onnx_probs.txt',"w") as f:
        for p in onnx:
            _p = float(p[1])
            f.write(str(_p) + '\n')
        f.close()

def main():
    
    if model_key == 'LightGBM': 
        train_LightGBM()
        convert_ONNX_LightGBM()
    elif model_key == 'TensorFlowNN': 
        train_TensorFlowNN()
        convert_ONNX_TensorFlowNN()

    if debug_mode_key['SHOW_ACCULACY'] == 1:
        showAccuracyPrecision()

    if debug_mode_key['SHOW_PLOT'] == 1:
        plt.show()

if __name__ == "__main__":

    global model_key
    global infile_key 
    global load_key 
    global calc_key 
    global object_key 
    global hyper_param_key 
    global df_key 
    global train_params_key 
    global output_name_key
    global debug_mode_key

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

    global model
    global onnx_model
    
    getExternalParams()    
    preprocessingData()

    main()
