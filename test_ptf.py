def estimar_textura_solo(teor_areia, teor_silte, teor_argila):
    # Verificando se os valores estão dentro do intervalo válido (0 a 100)
    if not (0 <= teor_areia <= 100) or not (0 <= teor_silte <= 100) or not (0 <= teor_argila <= 100):
        raise ValueError("Os teores devem estar no intervalo de 0 a 100.")
    
    if teor_areia > 85 and teor_silte < 15 and teor_argila < 10:
        return "Areia"
    elif teor_areia < 50 and teor_silte > 80 and teor_argila < 12:
        return "Silte"
    elif teor_areia < 45 and teor_silte < 40 and teor_argila > 40:
        return "Argila"
    elif 50 <= teor_areia <= 70 and 10 <= teor_silte <= 50 and teor_argila < 27:
        return "Franco Arenoso"
    elif 20 <= teor_areia <= 50 and 50 <= teor_silte <= 80 and teor_argila < 27:
        return "Franco Siltoso"
    elif teor_areia < 45 and teor_silte < 50 and 27 <= teor_argila <= 40:
        return "Franco Argiloso"
    elif 20 <= teor_areia <= 50 and 20 <= teor_silte <= 50 and 20 <= teor_argila <= 50:
        return "Franco"
    elif teor_areia < 45 and 50 <= teor_silte <= 80 and teor_argila > 40:
        return "Silte Argiloso"
    elif 50 <= teor_areia <= 70 and 20 <= teor_silte <= 50 and teor_argila > 27:
        return "Silte Argiloso Arenoso"
    else:
        return "Classificação não encontrada"
