# MolNetEnhancerMod
MolNetEnhancerMod is an alternative to the non-working MolNetEnhancer. It uses a modified approximation to retrieve SMILES classification, supporting both [NPClassifier](https://pubs.acs.org/doi/10.1021/acs.jnatprod.1c00399) and ClassyFire for robust molecular networking and detailed structural insights in metabolomics data.

## Select what you need
### [A. For Map chemical class information using NPClassifier to mass spectral molecular networks](https://github.com/froz9/MolNetEnhancerMod/tree/main?tab=readme-ov-file#a-for-map-chemical-class-information-using-npclassifier-to-mass-spectral-molecular-networks-1)

### [B. For Map chemical class information using Classyfire to mass spectral molecular networks](https://github.com/froz9/MolNetEnhancerMod/tree/main?tab=readme-ov-file#b-for-map-chemical-class-information-using-classyfire-to-mass-spectral-molecular-networks-1)

## A. For Map chemical class information using NPClassifier to mass spectral molecular networks
This script is a modification of the original [MolNetEnhancer](https://www.mdpi.com/2218-1989/9/7/144) workflow published by Madeleine Ernst. This version assigns chemical class annotations using [NPClassifier](https://pubs.acs.org/doi/10.1021/acs.jnatprod.1c00399).

### To map chemical class information to a mass spectral molecular network, you need to:

- Create a molecular network using the [classical](https://ccms-ucsd.github.io/GNPSDocumentation/quickstart/) or [feature-based](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/) workflow through the Global Natural Products Social Molecular Networking (GNPS) platform
- Perform in silico structure annotation using Network Annotation Propagation (NAP), DEREPLICATOR, or another tool of preference for in silico structure annotation

Then execute the code in [MolNenEnhancerMod_NPC.R](https://github.com/froz9/MolNetEnhancerMod/blob/main/MolNenEnhancerMod_NPC.R) line by line.
#### Disclaimer: Currently, the code only works with Feature Base Molecular Networking and NAP. DEREPLICATOR will be added soon.
 The only things that you need to replace are:

1. Your GNPS Job ID

![image](https://github.com/user-attachments/assets/85378ed6-fabc-405d-8d0f-473970502cb4)

2. Your NAP Job ID

![image](https://github.com/user-attachments/assets/c673eaca-ccc3-4ab2-8a89-902cc8efb775)

3. Replace the aforementioned ID's in the R script (To replicate results storage in [Madeleine GitHub](https://github.com/madeleineernst/RMolNetEnhancer/blob/master/Example_notebooks/ChemicalClasses_2_Network_FeatureBased.ipynb) I will use the same data)

4. Download the [MolNenEnhancerMod_NPC.R](https://github.com/froz9/MolNetEnhancerMod/blob/main/MolNenEnhancerMod_NPC.R)

![Captura de pantalla 2025-02-27 171412](https://github.com/user-attachments/assets/5e10169b-b663-4bd7-9ab0-e302032ed3d0)

5. Open it in your preferred R development environment, in this case, I use [RStudio](https://posit.co/downloads/) 

![Captura de pantalla 2025-02-27 171748](https://github.com/user-attachments/assets/07692410-39f2-4b0d-a484-84d21fee0319)

6. Replace the ID's from your work in the corresponding place

![Captura de pantalla 2025-02-27 172036](https://github.com/user-attachments/assets/c998e20e-b714-4eab-ae30-fe2968aa7c75)

7. Run the code as usual

8. After finished, it will generate a CSV file called **Molnetenhancer** in your working directory

![Captura de pantalla 2025-02-27 172333](https://github.com/user-attachments/assets/e7754806-2f9d-414d-bffa-369e52cc434c)

9. Open Cytoscape and import the Molecular Network File (You can find the graphML file inside of the folder called GNPS_output_graphML that was generated after run the previous code)
   ![Captura de pantalla 2025-02-27 172818](https://github.com/user-attachments/assets/c883467f-5192-4a13-ac63-65c250fd226d)

   The graphML file is located inside of **gnps_molecular_network_graphml**
   ![Captura de pantalla 2025-02-27 172847](https://github.com/user-attachments/assets/5f2e8c78-7ad8-4064-9d12-a7e0bd666b52)

   ![Captura de pantalla 2025-02-27 173349](https://github.com/user-attachments/assets/c22f97a2-0d9c-410a-b922-68b0419a0ba8)

   ![Captura de pantalla 2025-02-27 173428](https://github.com/user-attachments/assets/2b3753d2-4a0f-4660-9fad-982e70f25719)

10. After loaded, the graphML file will look like this:

     ![Captura de pantalla 2025-02-27 173805](https://github.com/user-attachments/assets/9780676b-f250-4388-9c58-1fd30d434ebb)

12. Import the **Molnethancer** CSV file into the network

    ![Captura de pantalla 2025-02-27 173945](https://github.com/user-attachments/assets/4296dfd0-8adf-44ea-acf3-fabfac93b0a4)

    ![Captura de pantalla 2025-02-27 194217](https://github.com/user-attachments/assets/5578bdb6-ceeb-4521-ac90-bfb9a3d9c931)

    **cluster.index** need to have the key

    ![Captura de pantalla 2025-02-27 194327](https://github.com/user-attachments/assets/a62ffc29-abb3-4b62-bbc2-93c50912fe9c)

    Then press **Ok**

14. Navigate to Style -> Fill Color

    ![Captura de pantalla 2025-02-27 194531](https://github.com/user-attachments/assets/82ee7554-e8b6-4f2f-9a83-98f082190a56)

    In **Column**, select the Classification level that you need (Pathway, Superclass, or Class)
    For the example, I chose **NPC_Pathway_Consensus**, in **Mapping Type** select **Discrete Mapping**
    The NPC Pathways will be displayed

    ![Captura de pantalla 2025-02-27 195019](https://github.com/user-attachments/assets/202a9917-dcae-4575-8ea4-0f3da3e809c3)

15. Select the colors that you want. For example, I right-clicked the in the row next to the **Alkaloids** and inside **Mapping Values Generator**, I selected **Rainbow**

    ![Captura de pantalla 2025-02-27 195433](https://github.com/user-attachments/assets/70808a75-c74c-4b78-899f-002cddf74746)

16. Your network enhanced with [NPClassifier](https://pubs.acs.org/doi/10.1021/acs.jnatprod.1c00399) is ready!

    ![Captura de pantalla 2025-02-27 195751](https://github.com/user-attachments/assets/d6f22299-4c3f-410e-a12e-acf43aa859ce)

## B. For Map chemical class information using Classyfire to mass spectral molecular networks
This script is a modification of the original [MolNetEnhancer](https://www.mdpi.com/2218-1989/9/7/144) workflow published by Madeleine Ernst. This version assigns chemical class annotations using [NPClassifier](https://pubs.acs.org/doi/10.1021/acs.jnatprod.1c00399).
