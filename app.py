from joblib import load
from os import listdir, path
from flask import Flask, render_template, request, send_file, jsonify
import matplotlib.pyplot as plt
import matplotlib

from numpy import array, mean
from io import BytesIO
from functions_ptf import PTF
from flask_cors import CORS

models = {}
resultados_globais = {}
files_in_folder = listdir("./models")
for model in files_in_folder:
    if model.endswith('.pkl'):
        name = model.replace('.pkl', '')
        path_file = path.join("./models", model)
        models[name] = load(path_file)

#The last, Agg, is a non-interactive backend that can only write to files
matplotlib.use('agg')
app = Flask(__name__)
CORS(app, origins='*')


def create_graph(resultados):
    ordem = sorted(resultados.keys())
    dicionario_ordenado = {chave: resultados[chave] for chave in ordem}

    resultados_cc = {k: v for k, v in dicionario_ordenado.items() if k.endswith('_cc')}
    resultados_pmp = {k: v for k, v in dicionario_ordenado.items() if k.endswith('_pmp')}

    modelos_cc = list(resultados_cc.keys())
    valores_cc = [valor.item() for valor in resultados_cc.values()]  # Converta os valores para escalares
    valores_pmp = [valor.item() for valor in resultados_pmp.values()]  # Converta os valores para escalares

    nomes = [name.replace('_cc', "") for name in modelos_cc]
    diferenca_valores = [(cc - pmp)*10 for cc, pmp in zip(valores_cc, valores_pmp)]

    # Crie um gráfico de barras com os resultados para _cc
    plt.figure(figsize=(9, 9))
    plt.bar(nomes, diferenca_valores, color='blue', label='AD')
    plt.axhline(y=mean(diferenca_valores), color='r', linestyle='--', linewidth=1, label=f'Média {round(mean(diferenca_valores), 2)}')
    for i, valor in enumerate(diferenca_valores):
        plt.text(i, valor + 0.01, str(round(valor, 2)), ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    plt.xlabel('Modelos')
    plt.ylabel('Resultados (mm/cm)')
    plt.title('Diferença entre Modelos CC e PMP')
    plt.xticks(rotation=90, ha='right')  # Ajuste a rotação e alinhamento dos rótulos
    plt.legend()  # Adicione a legenda
    plt.tight_layout()  # Ajuste o layout para evitar cortes

    # Salve o gráfico em um objeto buffer de bytes
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)
    plt.close()

    return buffer





def models_predict(inputs):
    results = {}
    for name, model in models.items():  # Iterar pelos pares (nome, modelo) no dicionário
        results[name] = float(model.predict(inputs))
    return results


def results(clay, silt, sand, density, porosity, org_mat):
    # Converta as strings em números float usando float()
    sand = float(sand)
    silt = float(silt)
    clay = float(clay)
    density = float(density)
    org_mat = float(org_mat)
    porosity = float(porosity)

    predict = [[clay, silt, sand, density, porosity, org_mat]]
    predict = array(predict)
    return predict


# @app.route('/', methods=['GET', 'POST'])
# def formulario():
#     global resultados_globais
#     if request.method == 'POST':
#         # Resgata os dados enviados pelo usuário
#         sand = request.form['sand']
#         clay = request.form['clay']
#         silt = request.form['silt']
#         density = request.form['density']
#         org_mat = request.form['org_mat']
#         porosity = request.form['porosity']
#         # Coloque na ordem clay - silt - sand - bulk_den - porosity - org_mat
#         resultados_globais = models_predict(results(clay, silt, sand, density, porosity, org_mat))
        
#         # textura = PTF.estimar_textura_solo(float(sand), float(silt), float(clay))
#         # return render_template('resultado.html', sand=sand, clay=clay, silt=silt, density=density, org_mat=org_mat, porosity=porosity, textura=textura)
#     return resultados_globais
#     # return render_template('formulario.html')

@app.route('/', methods=['GET'])
def home():
    return render_template('index.html')

@app.route('/ptf', methods=['POST'])
def send_ptfs():
    try:
        resultados_globais = {}
        data = request.get_json()
        print(data)
        required_fields = ['sand', 'clay', 'silt', 'density', 'org_mat', 'porosity']
        for field in required_fields:
            if field not in data:
                raise ValueError(f'O campo {field} é obrigatório.')
       
        if request.method == 'POST':
            # Resgata os dados enviados pelo usuário
            sand = data['sand']
            clay = data['clay']
            silt = data['silt']
            density = data['density']
            org_mat = data['org_mat']
            porosity = data['porosity']
            # Coloque na ordem clay - silt - sand - bulk_den - porosity - org_mat
            resultados_globais = models_predict(results(clay, silt, sand, density, porosity, org_mat))
            print(resultados_globais)
        else:
            resposta = {
                'status': 'error',
                'mensagem': 'Requisição mal formulada!'
            }
            return jsonify(resposta)
            
        resposta = {
            'status': 'success',
            'mensagem': 'Dados recebidos com sucesso',
            'dados': resultados_globais
        }
        
        return jsonify(resposta)
    
    except Exception as e:
        resposta = {
            'status': 'error',
            'mensagem': 'Erro ao processar os dados',
            'dados': str(e)
        }
        return jsonify(resposta),400
   



@app.route('/plot', methods=['GET', 'POST'])
def plot():
    global resultados_globais  # Acesse a variável global

    grafico_buffer = create_graph(resultados_globais)
    return send_file(grafico_buffer, mimetype='image/png')


if __name__ == '__main__':
    app.run(debug=True, threaded=False)
