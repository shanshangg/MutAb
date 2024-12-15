import argparse
import pandas as pd 
from autogluon.tabular import TabularDataset, TabularPredictor

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', default='input.csv')
    parser.add_argument('-o', '--output', default='output.csv')
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output

    df = pd.read_csv(input_file)
    mutant_info = df.iloc[:, 0]
    feat =  df.iloc[:, 1:97]
    feat = feat.rename(columns={0:'class'})

    predictor = TabularPredictor.load("..\\models\\TabularPredictor")

    pred_y = predictor.predict(feat)
    pred_prob_y = predictor.predict_proba(feat)
    
    output = pd.concat([mutant_info,pred_y,pred_prob_y],axis=1)
    output.to_csv(output_file,index=None)
    print("Done.")
