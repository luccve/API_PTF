import math
from typing import Optional

class PTF:
    @staticmethod
    def materiaOrganica(CO: float) -> float:
        return CO * 1.724
    
    @staticmethod
    def medrado(AT: float, SIL: float, ARG: float, CO: float, DS: float) -> list[float]:
        a1_2 = 0.13461391
        b1_2 = -1.58555722
        a2_2 = 0.04883605
        b2_2 = 0.19647777
        a3_2 = -0.00949548
        b3_2 = 0.56355803
        a4_2 = -0.00005212
        b4_2 = 1.43345189
        a5_2 = 0.01057297
        b5_2 = 1.14951659
        F2 = 0.01153532
        a2_1 = -0.01831805
        b2_1 = 0.89935543
        a3_1 = -0.01131157
        b3_1 = 1.00134021
        a4_1 = -0.00684340
        b4_1 = 1.11564515
        a5_1 = 0.01622120
        b5_1 = 0.48555009
        F1 = 2.01973831
        theta_p = 0.7022
        MO = PTF.materiaOrganica(CO)
        a1_4 = 1.82865464
        b1_4 = -0.57582713
        a2_4 = -0.04262312
        b2_4 = 0.59924802
        a3_4 = 1.51361637
        b3_4 = 0.00985830
        a4_4 = 0.67293635
        b4_4 = 0.05966974
        a5_4 = -0.00163216
        b5_4 = 2.66119588
        F4 = -0.64209861
        a1_3 = -0.00351317
        b1_3 = 2.42317427
        a2_3 = -0.03460735
        b2_3 = 0.58034639
        a3_3 = -0.01845649
        b3_3 = 0.85117517
        a4_3 = -0.04897867
        b4_3 = 0.59621217
        a5_3 = 0.00001279
        b5_3 = 5.38486557
        F3 = 1.58532681

        theta_s = (
            (a2_1 * ARG ** b2_1)
            + (a3_1 * AT ** b3_1)
            + (a4_1 * SIL ** b4_1)
            + (a5_1 * MO ** b5_1)
            + F1
        ) * theta_p

        theta_r = (
            a1_2 * DS ** b1_2
            + a2_2 * ARG ** b2_2
            + a3_2 * AT ** b3_2
            + a4_2 * SIL ** b4_2
            + a5_2 * MO ** b5_2
            + F2
        )

        n = (
            a1_4 * DS ** b1_4
            + a2_4 * ARG ** b2_4
            + a3_4 * AT ** b3_4
            + a4_4 * SIL ** b4_4
            + a5_4 * MO ** b5_4
            + F4
        )

        log_alpha = (
            a1_3 * DS ** b1_3
            + a2_3 * ARG ** b2_3
            + a3_3 * AT ** b3_3
            + a4_3 * SIL ** b4_3
            + a5_3 * MO ** b5_3
            + F3
        )

        return [log_alpha, n, theta_r, theta_s]
    
    @staticmethod
    def wenceslau(AT: float, SIL: float, ARG: float) -> float:
        eq = (
            1
            + 0.3591
            * (
                -0.02128887 * AT
                + -0.01005814 * SIL
                + -0.01901894 * ARG
                + 0.0001171219 * AT * SIL
                + 0.0002073924 * AT * ARG
                + 0.00006118707 * SIL * ARG
                + -0.000006373789 * AT * SIL * ARG
            )
        )
        AD = eq ** 2.784741
        return AD * 10

    @staticmethod
    def masutti(ARG: float, SIL: float) -> list[float]:
        cc_33 = (-1.5691 + 0.4289 * (ARG + SIL)) / 100
        pmp_1500 = (0.530482 + 0.301235 * SIL + 0.092822 * ARG) / 100
        return [cc_33, pmp_1500]
    
    @staticmethod
    def wenceslau(AT: float, SIL: float, ARG: float) -> float:
        eq = 1 + 0.3591 * (-0.02128887 * AT + -0.01005814 * SIL + -0.01901894 * ARG +
                         0.0001171219 * AT * SIL + 0.0002073924 * AT * ARG + 0.00006118707 * SIL * ARG +
                         -0.000006373789 * AT * SIL * ARG)
        AD = eq ** 2.784741
        return AD * 10

    @staticmethod
    def masutti(ARG: float, SIL: float) -> list[float]:
        cc_33 = (-1.5691 + 0.4289 * (ARG + SIL)) / 100
        pmp_1500 = (0.530482 + 0.301235 * SIL + 0.092822 * ARG) / 100
        return [cc_33, pmp_1500]

    @staticmethod
    def peraza(ARG: float, SIL: float, CO: float) -> list[float]:
        MO = PTF.materiaOrganica(CO) * 10

        cc_10 = 0.04947 + (ARG * 10) * 0.000377 + (SIL * 10) * 0.000343 + 0.001890 * MO
        cc_33 = 0.00767 + 0.000304 * (ARG * 10) + 0.000050 * (SIL * 10) + 0.001750 * MO
        pmp_1500 = -0.02310 + (ARG * 10) * 0.000298 + (SIL * 10) * 0.000176 + 0.000867 * MO
        return [cc_10, cc_33, pmp_1500]

    @staticmethod
    def reichert(ARG: float, SIL: float, CO: float, DS: float) -> list[float]:
        MO = PTF.materiaOrganica(CO) / 100
        cc_10 = 0.268 + 0.05 * (ARG / 100) + 0.24 * (SIL / 100 + ARG / 100) + 0.85 * MO - 0.127 * DS
        cc_33 = 0.106 + 0.29 * (ARG / 100 + SIL / 100) + 0.93 * MO - 0.048 * DS
        pmp_1500 = -0.04 + 0.15 * (ARG / 100) + 0.17 * (ARG / 100 + SIL / 100) + 0.91 * MO + 0.026 * DS
        return [cc_10, cc_33, pmp_1500]

    @staticmethod
    def oliveira(AT: float, SIL: float, ARG: float, DS: float) -> float:
        AD = (-0.000021 * AT * 10 + 0.000203 * SIL * 10 + 0.000054 * ARG * 10 + 0.021656 * DS) * DS
        return AD * 10

    @staticmethod
    def oliveira_p(AT: float, SIL: float, ARG: float, DS: Optional[float] = None) -> list[float]:
        cc_33 = (0.000333 * SIL * 10 + 0.000387 * ARG * 10) * DS
        pmp_1500 = (0.000038 * AT * 10 + 0.000153 * SIL * 10 + 0.000341 * ARG * 10 - 0.030861 * DS) * DS
        return [cc_33, pmp_1500]
    
    @staticmethod
    def estimar_textura_solo(teor_areia: float, teor_silte: float, teor_argila: float)-> str:
        # Verificando se os valores estão dentro do intervalo válido (0 a 100)
        if not (0 <= teor_areia <= 100) or not (0 <= teor_silte <= 100) or not (0 <= teor_argila <= 100):
            raise ValueError("Os teores devem estar no intervalo de 0 a 100.")
        
        if teor_areia > 70 and teor_silte < 15 or teor_argila < 10:
            return "Arenoso"
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

    @staticmethod
    def barros_simplificada(AT: float, SIL: float, ARG: float) -> list[float]:
        theta_s = 0.434714 + 0.114177 * (AT / 100) + 0.117845 * (SIL / 100)
        theta_r = 0.128617 + 0.14836 * (AT / 100) + 0.35705 * (ARG / 100)
        log_alpha = -1.07329 - 1.9578 * (ARG / 100)
        alpha = 10 ** log_alpha
        n = 1.134153 + 0.722216 * (AT / 100) + 0.39574 * (SIL / 100)

        return [alpha, n, theta_r, theta_s]

    @staticmethod
    def van_genuchten(alpha: float, n: float, theta_r: float, theta_s: float) -> float:
        def theta(alpha: float, n: float, theta_r: float, theta_s: float, potential: float) -> float:
            return theta_r + (theta_s - theta_r) / (1 + abs(alpha * potential) ** n) ** (1 - 1 / n)

        AD10 = theta(alpha, n, theta_r, theta_s, 10)
        AD33 = theta(alpha, n, theta_r, theta_s, 33)
        AD1500 = theta(alpha, n, theta_r, theta_s, 1500)

        return AD10 - AD1500

    @staticmethod
    def van_genuchtenList(alpha: float, n: float, theta_r: float, theta_s: float) -> list[dict[str, float]]:
        def theta(alpha: float, n: float, theta_r: float, theta_s: float, potential: float) -> float:
            return theta_r + (theta_s - theta_r) / (1 + abs(alpha * potential) ** n) ** (1 - 1 / n)

        list_result = []
        potenciais = [1, 6, 10, 33, 100, 500, 1000, 1100, 1200, 1300, 1400, 1500]

        for i in potenciais:
            list_result.append({"index": i, "theta": theta(alpha, n, theta_r, theta_s, i)})

        return list_result
    
    @staticmethod
    def wenceslau(AT: float, SIL: float, ARG: float) -> float:
        eq = (1 + 0.3591 * ((-0.02128887 * AT) + (-0.01005814 * SIL) + (-0.01901894 * ARG) + (0.0001171219 * AT * SIL) + (0.0002073924 * AT * ARG) + (0.00006118707 * SIL * ARG) + (-0.000006373789 * AT * SIL * ARG)))
        AD = eq ** 2.784741
        return AD * 10

    @staticmethod
    def masutti(ARG: float, SIL: float) -> list[float]:
        cc_33 = ((-1.5691 + 0.4289 * (ARG + SIL)) / 100)
        pmp_1500 = ((0.530482 + 0.301235 * SIL + 0.092822 * ARG) / 100)
        return [cc_33, pmp_1500]

    @staticmethod
    def peraza(ARG: float, SIL: float, CO: float) -> list[float]:
        MO = PTF.materiaOrganica(CO) * 10

        cc_10 = 0.04947 + (ARG * 10) * 0.000377 + (SIL * 10) * 0.000343 + 0.001890 * MO
        cc_33 = 0.00767 + 0.000304 * (ARG * 10) + 0.000050 * (SIL * 10) + 0.001750 * MO
        pmp_1500 = -0.02310 + (ARG * 10) * 0.000298 + (SIL * 10) * 0.000176 + 0.000867 * MO
        return [cc_10, cc_33, pmp_1500]

    @staticmethod
    def reichert(ARG: float, SIL: float, CO: float, DS: float) -> list[float]:
        MO = PTF.materiaOrganica(CO) / 100
        cc_10 = 0.268 + 0.05 * (ARG / 100) + 0.24 * (SIL / 100 + ARG / 100) + 0.85 * MO - 0.127 * DS
        cc_33 = 0.106 + 0.29 * (ARG / 100 + SIL / 100) + 0.93 * MO - 0.048 * DS
        pmp_1500 = -0.04 + 0.15 * (ARG / 100) + 0.17 * (ARG / 100 + SIL / 100) + 0.91 * MO + 0.026 * DS
        return [cc_10, cc_33, pmp_1500]

    @staticmethod
    def oliveira(AT: float, SIL: float, ARG: float, DS: float) -> float:
        AD = (-0.000021 * AT * 10 + 0.000203 * SIL * 10 + 0.000054 * ARG * 10 + 0.021656 * DS) * DS
        return AD * 10

    @staticmethod
    def oliveira_p(AT: float, SIL: float, ARG: float, DS: float) -> list[float]:
        cc_33 = (0.000333 * SIL * 10 + 0.000387 * ARG * 10) * DS
        pmp_1500 = (0.000038 * AT * 10 + 0.000153 * SIL * 10 + 0.000341 * ARG * 10 - 0.030861 * DS) * DS
        return [cc_33, pmp_1500]

    @staticmethod
    def vandenBerg(SIL: float, ARG: float, MO: float) -> list[float]:
        carbonoOrganico = MO * 0.58
        CC_10 = (10.88 + 0.347 * ARG + 0.211 * SIL + 0.1756 * carbonoOrganico) / 100
        PMP = (3.83 + 0.272 * ARG + 0.212 * SIL) / 100
        return [CC_10, PMP]

    @staticmethod
    def vanderBerg_vanGenuchten(ARG: float, AT: float, SIL: float, CO: float, CTC: float) -> list[float]:
        n = 1 / (1 - 0.503 - 0.0027 * (SIL + ARG) - 0.066 * CO + 0.0094 + CTC)
        theta_s = (84.1 - 0.206 * ARG - 0.322 * (AT + SIL)) / 100
        theta_r = (0.308 * ARG) / 100
        log_alpha = -0.627
        alpha = 10 ** log_alpha
        return [alpha, n, theta_r, theta_s]

    @staticmethod
    def andrade_stone(SIL: float, ARG: float, AT: float, CO: float) -> list[float]:
        MO = PTF.materiaOrganica(CO)
        U_cc_6 = 0.0019 * ARG + 0.0024 * SIL + 0.2143
        return [U_cc_6, 0]

    @staticmethod
    def hodnett(ARG: float, AT: float, SIL: float, CTC: float, pH: float, CO: float, DS: float) -> list[float]:
        theta_s = (82.072 + (0.089 * ARG) - (31.357 * DS) + (0.027 * CTC) + (0.517 * pH) - (0.0006 * AT * ARG)) / 100
        theta_r = (23.133 + (-0.172 * AT) + (0.211 * CTC) + (-0.849 * pH) + (0.0012 * (ARG * ARG)) + (0.0029 * AT * ARG)) / 100
        log_alpha = (-4.237 - (3.423 * SIL) + (4.288 * CO) - (0.801 * CTC) - (11.07 * pH) + (0.027 * SIL * SIL)) / 100
        n = (67.093 + (-0.907 * ARG) + (-0.574 * CO) + (1.396 * pH) + (0.0056 * (ARG * ARG))) / 100
        alpha = 10 ** log_alpha
        return [alpha, n, theta_r, theta_s]
    
    @staticmethod
    def giarola(SIL: float, ARG: float) -> list[float]:
        cc_10 = 0.081 + 0.005 * SIL + 0.004 * ARG
        pmp_1500 = -0.031 + 0.005 * SIL + 0.003 * ARG
        return [cc_10, pmp_1500]

    @staticmethod
    def DASAM(AT: float, SIL: float, ARG: float, DS: float, CO: float) -> list[float]:
        MO = 16 / 1000
        log_alpha = 0.8118 + (0.8861 * AT / 100) + (-1.1907 * ARG / 100) + (-1.5140 * DS)

        alpha = 10 ** log_alpha
        n = 1.1527 + (-5.5341 * MO) + (0.7427 * AT / 100) + (0.4135 * SIL / 100)
        theta_r = 0.0858 + (1.1846 * MO) + (-0.1671 * AT / 100) + (0.3516 * ARG / 100) + (0.0290 * DS)
        theta_s = 1 + (-0.370000 * DS)
        return [alpha, n, theta_r, theta_s]
    
    @staticmethod
    def giarolab(SIL: float, ARG: float) -> list[float]:
        cc_10 = 0.081 + 0.005 * SIL + 0.004 * ARG
        pmp_1500 = 0.024 + 0.005 * SIL + 0.003 * ARG
        return [cc_10, pmp_1500]

    @staticmethod
    def densidade(AT: float, ARG: float, CARB: float) -> float:
        bd = 1.286 + 3.208 * 0.001 * AT - 2.013 * 0.001 * ARG
        bd_oc = 1.358 + 2.79 * 0.001 * AT - 2.328 * 0.001 * ARG - 0.0052 * CARB
        return bd_oc if CARB is not None else bd

    @staticmethod
    def tomasella_n1(
        AG: float,
        AF: float,
        SIL: float,
        ARG: float,
        DS: float,
        EQ: float
    ) -> list[float]:
        log_alpha = (
            264.4658
            + 1.2123 * ARG
            - 378.6112 * EQ
            - 328.3456 * DS
            + 0.0052 * AG * AF
            + 0.00775 * AG * SIL
            + 0.0963 * AF * ARG
            + 0.0616 * AG * AG
        )
        n = (
            2.1909
            - 1.5296 * EQ
            - 0.000299 * AG * SIL
            - 0.000345 * AF * ARG
            - 0.000105 * AG * AG
            + 0.000025 * AF * AF
        )

        alpha = 10 ** (log_alpha)

        theta_s = (
            0.8219
            - 0.000177 * SIL
            + 2.2324 * EQ
            - 0.2867 * DS
            + 0.000049 * AG * SIL
            - 0.000029 * AG * ARG
            + 0.000027 * AF * ARG
            - 0.000008 * AG * AG
        )

        theta_r = (
            -0.1336
            + 0.0025 * SIL
            + 0.0034 * ARG
            + 0.3991 * EQ
            + 0.0768 * DS
            - 0.000048 * SIL * SIL
            - 0.000013 * ARG * ARG
        )

        return [alpha, n, theta_r, theta_s]

    @staticmethod
    def tomasella_n2(
        AG: float,
        AF: float,
        SIL: float,
        ARG: float,
        CO: float,
        UEQ: float,
    ) -> list[float]:
        carbonoOrganico = CO  # % to g/kg
        theta_s = (
            42.1253
            + (57.9073 * UEQ)
            + (1.5587 * carbonoOrganico)
            + (-0.0039 * AG * AF)
            + (0.0033 * AG * SIL)
            + (-0.0041 * AG * ARG)
            + (-0.0042 * SIL * ARG)
        ) / 100

        theta_r = (
            -1.0223
            + (0.1754 * SIL)
            + (0.2770 * ARG)
            + (48.3545 * UEQ)
            + (-1.5731 * carbonoOrganico)
            + (0.0014 * AF * ARG)
            + (-0.0046 * SIL ** 2)
            + (-0.0013 * ARG ** 2)
        ) / 100

        log_alpha = (
            -200.8303
            + (1.0190 * ARG)
            + (184.9474 * UEQ)
            + (-0.0284 * AG * AF)
            + (0.0378 * AG * SIL)
            + (0.0634 * AF * ARG)
            + (-0.0625 * SIL * ARG)
            + (0.0657 * AG ** 2)
        )

        alpha = 10 ** (log_alpha / 100) * 0.75
        n = (
            217.4984
            + (-147.4518 * UEQ)
            + (-0.0295 * AG * SIL)
            + (-0.0375 * AF * ARG)
            + (-0.0099 * AG ** 2)
            + (0.007 * AF ** 2)
        ) / 100
        return [alpha, n, theta_r, theta_s]
    
    @staticmethod
    def tomasella_n3(
        AG: float,
        AF: float,
        SIL: float,
        ARG: float,
        CO: float,
        DS: float,
    ) -> list[float]:
        carbonoOrganico = CO * 10  # % to g/kg
        theta_s = (
            91.6203
            + (-30.0046 * DS)
            + (1.5925 * carbonoOrganico)
            + (0.0022 * AG * SIL)
            + (-0.0036 * AG * ARG)
            + (-0.0018 * AG ** 2 + (-0.001 * AF ** 2))
        ) / 100

        theta_r = (
            23.3867
            + (0.1103 * ARG)
            + (-4.7949 * DS)
            + (0.0047 * SIL * ARG)
            + (-0.0027 * AG ** 2)
            + (-0.0022 * AF ** 2)
            + (-0.0048 * SIL ** 2)
        ) / 100

        log_alpha = (
            205.6546
            + (-2.556 * SIL)
            + (-0.1329 * ARG)
            + (-247.4904 * DS)
            + (-0.0189 * AG * AF)
            + (0.1177 * AG * SIL)
            + (0.0517 * AF * ARG)
            + (0.0617 * (AG ** 2))
        )

        alpha = 10 ** (log_alpha / 100) * 0.75
        n = ((168.8617 + (-0.0258 * AG * SIL) + (-0.0261 * AF * ARG) + (0.0093 * AF ** 2) + (-0.0077 * SIL ** 2))) / 100
        return [alpha, n, theta_r, theta_s]

    @staticmethod
    def tomasella_n4(
        AG: float,
        AF: float,
        SIL: float,
        ARG: float,
        CO: float,
    ) -> list[float]:
        theta_s = (
            36.9
            + (0.37 * SIL)
            + (3.26 * CO)
            + (-0.002 * AG * ARG)
            + (0.003 * AF * ARG)
            + (-0.003 * SIL * ARG)
            + (0.003 * ARG * ARG)
        ) / 100

        theta_r = (
            (15.76 + (0.14 * ARG) + (0.005 * SIL * ARG) + (-0.003 * (AG ** 2)) + (-0.002 * AF ** 2))
            + (-0.005 * SIL ** 2)
        ) / 100

        log_alpha = (-237.01 + (3.62 * AG) + (0.004 * AG * SIL) + (0.09 * AF * ARG) + (0.018 * (ARG ** 2))) / 100
        alpha = 10 ** log_alpha
        n = (170.63 + (-0.018 * AG * SIL) + (-0.031 * AF * ARG) + (0.009 * AF ** 2 + (-0.008 * SIL ** 2))) / 100
        return [alpha, n, theta_r, theta_s]

    @staticmethod
    def leonorAssad(AT: float) -> float:
        AD = (12.76278562 - 0.0000098726 * AT ** 3) / 100
        return AD * 10

    @staticmethod
    def arruda(SIL: float, ARG: float) -> list[float]:
        cc_33 = (
            3.07439 + 0.629239 * (SIL + ARG) - 0.00343813 * (SIL + ARG) ** 2
        ) / 100
        pmp = (398.889 * (SIL + ARG)) / (1308.09 + (SIL + ARG)) / 100
        return [cc_33, pmp]