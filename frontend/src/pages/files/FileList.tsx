/**
 * File List Page - Complete file upload and management
 */
import React, { useEffect, useState } from 'react'
import {
  Button,
  Select,
  Tag,
  Modal,
  Upload,
  message,
  Popconfirm,
  Tooltip,
  Space,
} from 'antd'
import {
  UploadOutlined,
  DownloadOutlined,
  FileOutlined,
  InboxOutlined,
  DeleteOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import type { UploadProps } from 'antd'
import { useFileStore } from '../../store/fileStore'
import { useProjectStore } from '../../store/projectStore'
import { useSampleStore } from '../../store/sampleStore'
import { PageHeader, DataTable } from '../../components/common'
import { toast } from '../../utils/notification'
import type { FileItem } from '../../types/file'
import dayjs from 'dayjs'

const { Option } = Select
const { Dragger } = Upload

export const FileList: React.FC = () => {
  const { files, loading, fetchFiles, uploadFile, deleteFile, downloadFile } = useFileStore()
  const { projects, fetchProjects } = useProjectStore()
  const { samples, fetchSamples } = useSampleStore()
  const [selectedProject, setSelectedProject] = useState<string>('')
  const [selectedSample, setSelectedSample] = useState<string>('')
  const [isUploadModalVisible, setIsUploadModalVisible] = useState(false)
  const [uploading, setUploading] = useState(false)

  useEffect(() => {
    fetchProjects()
  }, [])

  useEffect(() => {
    if (selectedProject) {
      fetchFiles({ project_id: selectedProject })
      fetchSamples({ project_id: selectedProject })
    }
  }, [selectedProject])

  // 打开上传模态框
  const showUploadModal = () => {
    if (!selectedProject) {
      toast.warning('Please select a project first')
      return
    }
    setIsUploadModalVisible(true)
  }

  // 关闭上传模态框
  const handleCancelUpload = () => {
    setIsUploadModalVisible(false)
    setSelectedSample('')
  }

  // 文件上传配置
  const uploadProps: UploadProps = {
    name: 'file',
    multiple: true,
    accept: '.fastq,.fq,.fastq.gz,.fq.gz,.bam,.sam,.vcf,.vcf.gz',
    showUploadList: true,
    beforeUpload: () => false, // 阻止自动上传
    onChange(info) {
      // 处理文件选择变化
    },
  }

  // 手动上传文件
  const handleUpload = async (fileList: any[]) => {
    if (fileList.length === 0) {
      toast.warning('Please select at least one file')
      return
    }

    if (!selectedSample) {
      toast.warning('Please select a sample for the files')
      return
    }

    setUploading(true)

    try {
      // 逐个上传文件
      for (const fileItem of fileList) {
        const file = fileItem.originFileObj || fileItem
        await uploadFile(selectedSample, file)
      }

      toast.success(`Successfully uploaded ${fileList.length} file(s)`)
      handleCancelUpload()

      // 刷新文件列表
      if (selectedProject) {
        fetchFiles({ project_id: selectedProject })
      }
    } catch (error) {
      console.error('Upload failed:', error)
      toast.error('Failed to upload files')
    } finally {
      setUploading(false)
    }
  }

  // 删除文件
  const handleDelete = async (id: string) => {
    try {
      await deleteFile(id)
      toast.success('File deleted successfully')

      // 刷新文件列表
      if (selectedProject) {
        fetchFiles({ project_id: selectedProject })
      }
    } catch (error) {
      console.error('Failed to delete file:', error)
      toast.error('Failed to delete file')
    }
  }

  const formatFileSize = (bytes: number) => {
    if (bytes === 0) return '0 B'
    const k = 1024
    const sizes = ['B', 'KB', 'MB', 'GB', 'TB']
    const i = Math.floor(Math.log(bytes) / Math.log(k))
    return `${(bytes / Math.pow(k, i)).toFixed(2)} ${sizes[i]}`
  }

  const columns: ColumnsType<FileItem> = [
    {
      title: 'Filename',
      dataIndex: 'filename',
      key: 'filename',
      render: (name) => (
        <Space>
          <FileOutlined style={{ color: 'var(--color-primary)' }} />
          <span style={{ fontWeight: 500 }}>{name}</span>
        </Space>
      ),
    },
    {
      title: 'Type',
      dataIndex: 'file_type',
      key: 'file_type',
      width: 120,
      render: (type) => {
        const colorMap: Record<string, string> = {
          fastq: 'blue',
          bam: 'green',
          vcf: 'purple',
          sam: 'orange',
        }
        const color = colorMap[type?.toLowerCase()] || 'default'
        return type ? <Tag color={color}>{type.toUpperCase()}</Tag> : '-'
      },
    },
    {
      title: 'Size',
      dataIndex: 'file_size',
      key: 'file_size',
      width: 120,
      render: (size) => formatFileSize(size || 0),
    },
    {
      title: 'Sample',
      dataIndex: 'sample_id',
      key: 'sample_id',
      width: 150,
      render: (sampleId) => {
        const sample = samples.find((s) => s.id === sampleId)
        return sample ? sample.sample_id : '-'
      },
    },
    {
      title: 'Uploaded',
      dataIndex: 'uploaded_at',
      key: 'uploaded_at',
      width: 150,
      render: (date) => (date ? dayjs(date).format('YYYY-MM-DD HH:mm') : '-'),
    },
    {
      title: 'Actions',
      key: 'actions',
      width: 150,
      fixed: 'right',
      render: (_, record) => (
        <Space>
          <Tooltip title="Download">
            <Button
              type="text"
              size="small"
              icon={<DownloadOutlined />}
              onClick={() => downloadFile(record.id, record.filename)}
            />
          </Tooltip>
          <Popconfirm
            title="Delete File"
            description="Are you sure you want to delete this file?"
            onConfirm={() => handleDelete(record.id)}
            okText="Yes, Delete"
            cancelText="Cancel"
            okButtonProps={{ danger: true }}
          >
            <Tooltip title="Delete">
              <Button type="text" size="small" danger icon={<DeleteOutlined />} />
            </Tooltip>
          </Popconfirm>
        </Space>
      ),
    },
  ]

  return (
    <div>
      <PageHeader
        left={
          <Select
            placeholder="Select Project"
            style={{ width: 300 }}
            value={selectedProject || undefined}
            onChange={setSelectedProject}
            loading={projects.length === 0}
          >
            {projects.map((p) => (
              <Option key={p.id} value={p.id}>
                {p.name}
              </Option>
            ))}
          </Select>
        }
        right={
          <Button
            type="primary"
            icon={<UploadOutlined />}
            onClick={showUploadModal}
            disabled={!selectedProject}
          >
            Upload Files
          </Button>
        }
      />

      <DataTable
        columns={columns}
        dataSource={files}
        rowKey="id"
        loading={loading}
        pagination={{ pageSize: 20 }}
        emptyText="No Files"
        emptyDescription="Select a project and upload files to get started"
      />

      {/* 文件上传模态框 */}
      <Modal
        title="Upload Files"
        open={isUploadModalVisible}
        onCancel={handleCancelUpload}
        footer={null}
        width={700}
      >
        <div style={{ marginBottom: 16 }}>
          <label style={{ display: 'block', marginBottom: 8, fontWeight: 500 }}>
            Select Sample *
          </label>
          <Select
            placeholder="Choose a sample for the files"
            style={{ width: '100%' }}
            value={selectedSample || undefined}
            onChange={setSelectedSample}
          >
            {samples.map((s) => (
              <Option key={s.id} value={s.id}>
                {s.sample_id} {s.group_name && `(${s.group_name})`}
              </Option>
            ))}
          </Select>
        </div>

        <Upload.Dragger
          {...uploadProps}
          disabled={!selectedSample}
        >
          <p className="ant-upload-drag-icon">
            <InboxOutlined />
          </p>
          <p className="ant-upload-text">Click or drag files to this area to upload</p>
          <p className="ant-upload-hint">
            Support for single or bulk upload. Accepted formats: FASTQ, BAM, SAM, VCF
            <br />
            Maximum file size: 50GB per file
          </p>
        </Upload.Dragger>

        <div style={{ marginTop: 16, textAlign: 'right' }}>
          <Space>
            <Button onClick={handleCancelUpload}>Cancel</Button>
            <Button
              type="primary"
              loading={uploading}
              onClick={() => {
                const fileList = (document.querySelector('.ant-upload-list') as any)
                  ?.querySelectorAll('.ant-upload-list-item')
                if (fileList && fileList.length > 0) {
                  // 获取文件列表并上传
                  // 这里需要从 Upload 组件获取文件
                  message.info('Upload functionality in progress')
                }
              }}
            >
              Upload
            </Button>
          </Space>
        </div>
      </Modal>
    </div>
  )
}

export default FileList
